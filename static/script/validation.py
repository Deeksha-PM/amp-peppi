"""Sequence cleaning, FASTA parsing and input validation (no Streamlit dependency).

Used by both the web app (pages/1_Prediction.py) and the command-line tool (predict.py),
so the length limits and residue checks are enforced identically in both.
"""

import re
from dataclasses import dataclass, field

import pandas as pd

from feature_extraction import (
    AMINO_ACIDS,
    AMP_MAX_LEN,
    AMP_MIN_LEN,
    PROTEIN_MAX_LEN,
    PROTEIN_MIN_LEN,
)

LIMITS = {
    "AMP": (AMP_MIN_LEN, AMP_MAX_LEN),
    "Protein": (PROTEIN_MIN_LEN, PROTEIN_MAX_LEN),
}

# Upper bound on all-vs-all combinations, to keep a single request from running for ever.
MAX_BULK_PAIRS = 5000


def clean_sequence(seq) -> str:
    """Upper-case the sequence and drop whitespace and any non-letter characters.

    Invalid amino-acid letters are deliberately kept so that validation can report them.
    """
    seq = str(seq).upper().strip().replace(" ", "")
    return re.sub(r"[^A-Z]", "", seq)


def check_sequence(raw, label: str):
    """Clean one sequence and check it against the rules for `label` ("AMP" or "Protein").

    Returns (cleaned sequence, list of problems, list of removed non-letter characters).
    An empty problem list means the sequence is valid.
    """
    min_len, max_len = LIMITS[label]
    if raw is None or pd.isna(raw):
        raw = ""
    compact = re.sub(r"\s+", "", str(raw))
    cleaned = clean_sequence(raw)
    removed = sorted(set(re.findall(r"[^A-Za-z]", compact)))

    if not cleaned:
        return cleaned, ["the sequence is empty."], removed

    problems = []
    invalid = sorted(set(cleaned) - set(AMINO_ACIDS))
    if invalid:
        problems.append(
            f"contains unsupported residue(s) {', '.join(invalid)}. "
            f"Only the 20 standard amino acids ({AMINO_ACIDS}) are allowed."
        )
    if not (min_len <= len(cleaned) <= max_len):
        problems.append(
            f"length is {len(cleaned)}, which is outside the allowed range of "
            f"{min_len}-{max_len} residues."
        )
    return cleaned, problems, removed


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Return [(id, raw sequence), ...] from FASTA text; raise ValueError if it is not FASTA."""
    text = str(text).lstrip("﻿")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("The input is empty.")
    if not lines[0].startswith(">"):
        raise ValueError(
            "This does not look like FASTA: the first non-empty line must be a header "
            "starting with '>'."
        )

    records, current_id, parts = [], None, []
    for line in lines:
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, "".join(parts)))
            header = line[1:].strip()
            current_id = header.split()[0] if header else f"seq_{len(records) + 1}"
            parts = []
        else:
            parts.append(line)
    records.append((current_id, "".join(parts)))
    return records


def count_fasta_records(text: str) -> int:
    return sum(1 for line in str(text).splitlines() if line.lstrip().startswith(">"))


def parse_single_input(text: str, default_id: str):
    """Accept a plain sequence or pasted FASTA. Returns (id, raw sequence, warnings)."""
    text = str(text)
    if not text.lstrip().startswith(">"):
        return default_id, text, []
    records = parse_fasta(text)
    warnings = []
    if len(records) > 1:
        warnings.append(
            f"{default_id}: {len(records)} FASTA records were pasted; only the first "
            f"('{records[0][0]}') is used."
        )
    return records[0][0], records[0][1], warnings


def length_status(text: str, label: str) -> str:
    """One-line summary for the live length counter under a text box."""
    min_len, max_len = LIMITS[label]
    try:
        _, raw, _ = parse_single_input(text, label)
    except ValueError:
        return f"{label}: could not read the pasted FASTA."
    cleaned = clean_sequence(raw)
    n = len(cleaned)
    if n and not (min_len <= n <= max_len):
        msg = f"{n} residues - outside the allowed range of {min_len}-{max_len}"
    else:
        msg = f"{n} residues (allowed {min_len}-{max_len})"
    invalid = sorted(set(cleaned) - set(AMINO_ACIDS))
    if invalid:
        msg += f"; unsupported residue(s): {', '.join(invalid)}"
    return msg


@dataclass
class Validated:
    valid: list = field(default_factory=list)      # [{"id", "seq"}] cleaned, passed all checks
    errors: list = field(default_factory=list)     # human-readable messages
    invalid: list = field(default_factory=list)    # rows for the "invalid records" download
    warnings: list = field(default_factory=list)   # non-fatal notes


def validate_records(records, label: str) -> Validated:
    """Validate [(id, raw sequence), ...] for `label` ("AMP" or "Protein")."""
    result = Validated()
    if not records:
        result.errors.append(f"No {label} sequences were found in the input.")
        return result

    used_ids, dup_counts = set(), {}

    def unique_id(seq_id):
        if seq_id not in used_ids:
            used_ids.add(seq_id)
            return seq_id
        n = dup_counts.get(seq_id, 1)
        while f"{seq_id}_{n + 1}" in used_ids:
            n += 1
        n += 1
        dup_counts[seq_id] = n
        new_id = f"{seq_id}_{n}"
        used_ids.add(new_id)
        result.warnings.append(f"{label} ID '{seq_id}' is duplicated; renamed to '{new_id}'.")
        return new_id

    for seq_id, raw in records:
        seq_id = unique_id(seq_id)
        cleaned, problems, removed = check_sequence(raw, label)
        if removed and cleaned:
            result.warnings.append(
                f"{label} '{seq_id}': removed non-letter character(s) {' '.join(removed)}."
            )
        if problems:
            reason = " ".join(problems)
            result.errors.append(f"{label} '{seq_id}': {reason}")
            result.invalid.append(
                {"type": label, "id": seq_id, "length": len(cleaned),
                 "reason": reason, "sequence": cleaned}
            )
        else:
            result.valid.append({"id": seq_id, "seq": cleaned})
    return result


def validate_frame(
    frame: pd.DataFrame,
    amp_col: str = "amp_sequence",
    protein_col: str = "protein_sequence",
    max_errors: int = 20,
) -> pd.DataFrame:
    """Validate a CSV-style frame; return a copy with cleaned sequences or raise ValueError."""
    missing = [c for c in (amp_col, protein_col) if c not in frame.columns]
    if missing:
        raise ValueError(f"Input is missing required sequence columns: {missing}")

    cleaned_frame, errors = frame.copy(), []
    for column, label in ((amp_col, "AMP"), (protein_col, "Protein")):
        cleaned_values = []
        for row, raw in enumerate(frame[column], start=2):
            cleaned, problems, _ = check_sequence(raw, label)
            if problems:
                errors.append(f"{label} at input row {row}: {' '.join(problems)}")
            cleaned_values.append(cleaned)
        cleaned_frame[column] = cleaned_values

    if errors:
        shown = "\n".join(errors[:max_errors])
        extra = f"\n... and {len(errors) - max_errors} more." if len(errors) > max_errors else ""
        raise ValueError(f"Input validation failed:\n{shown}{extra}")
    return cleaned_frame
