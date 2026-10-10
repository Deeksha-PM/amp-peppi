"""Calculate the paired-sequence features expected by the saved model."""

import pickle
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import protpy


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Allowed sequence lengths (inclusive). Shared by the web app and the command-line tool.
AMP_MIN_LEN, AMP_MAX_LEN = 5, 50
PROTEIN_MIN_LEN, PROTEIN_MAX_LEN = 51, 1500

PHYS_FUNCTIONS = (
    "boman_index",
    "hydrophobic_moment",
    "aliphatic_index",
    "aromaticity",
    "gravy",
    "instability_index",
    "isoelectric_point",
    "molecular_weight",
    "shannon_entropy",
    "aggregation_propensity",
    "charge_distribution",
    "hydrophobic_polar_charged_composition",
)


def _normalise_sequence(sequence: object, column: str, row: int) -> str:
    if pd.isna(sequence):
        raise ValueError(f"Missing {column} at input row {row}.")
    normalised = "".join(str(sequence).split()).upper()
    if not normalised:
        raise ValueError(f"Empty {column} at input row {row}.")
    invalid = sorted(set(normalised) - set(AMINO_ACIDS))
    if invalid:
        raise ValueError(
            f"Unsupported residue(s) in {column} at input row {row}: {invalid}. "
            "Only the 20 standard amino-acid letters are supported."
        )
    return normalised


# Descriptors depend only on the sequence, so repeated sequences (e.g. the same protein
# paired with many AMPs, or the same input re-run) are computed once. The returned dict
# is shared between calls and must be treated as read-only.
@lru_cache(maxsize=4096)
def _describe_sequence(sequence: str) -> dict[str, float]:
    """Calculate AAC, CTD and PHYS values using the training feature definitions."""
    length = len(sequence)
    features = {f"AAC_{aa}": 100.0 * sequence.count(aa) / length for aa in AMINO_ACIDS}

    def add_block(result: pd.DataFrame, prefix: str) -> None:
        features.update({f"{prefix}_{name}": float(value)
                         for name, value in result.iloc[0].items()})

    add_block(protpy.ctd_(sequence, all_ctd=True), "CTD")
    for function_name in PHYS_FUNCTIONS:
        add_block(getattr(protpy, function_name)(sequence), "PHYS")

    if not np.isfinite(list(features.values())).all():
        raise ValueError("Descriptor calculation produced a non-finite feature value.")
    return features


def extract_features(
    frame: pd.DataFrame,
    model: object,
    amp_col: str = "amp_sequence",
    protein_col: str = "protein_sequence",
) -> pd.DataFrame:
    """Extract and order features from a CSV-loaded frame for this fitted pipeline."""
    missing_columns = [c for c in (amp_col, protein_col) if c not in frame.columns]
    if missing_columns:
        raise ValueError(f"Input is missing required sequence columns: {missing_columns}")
    if not hasattr(model, "feature_names_in_"):
        raise TypeError("The fitted model must expose feature_names_in_ for safe alignment.")

    required = list(model.feature_names_in_)
    amp_sequences = [
        _normalise_sequence(sequence, amp_col, row)
        for row, sequence in enumerate(frame[amp_col], start=2)
    ]
    protein_sequences = [
        _normalise_sequence(sequence, protein_col, row)
        for row, sequence in enumerate(frame[protein_col], start=2)
    ]
    unique_sequences = dict.fromkeys(amp_sequences + protein_sequences)
    descriptors = {sequence: _describe_sequence(sequence) for sequence in unique_sequences}
    values: dict[str, np.ndarray] = {}

    for feature in required:
        if feature.startswith("pep_"):
            descriptor = feature[len("pep_"):]
            values[feature] = np.asarray(
                [descriptors[sequence][descriptor] for sequence in amp_sequences]
            )
        elif feature.startswith("pro_"):
            descriptor = feature[len("pro_"):]
            values[feature] = np.asarray(
                [descriptors[sequence][descriptor] for sequence in protein_sequences]
            )
        elif feature.startswith("x_prod_") or feature.startswith("x_absdiff_"):
            prefix = "x_prod_" if feature.startswith("x_prod_") else "x_absdiff_"
            descriptor = feature[len(prefix):]
            pep_values = np.asarray(
                [descriptors[sequence][descriptor] for sequence in amp_sequences],
                dtype=float,
            )
            pro_values = np.asarray(
                [descriptors[sequence][descriptor] for sequence in protein_sequences],
                dtype=float,
            )
            values[feature] = (
                pep_values * pro_values
                if prefix == "x_prod_"
                else np.abs(pep_values - pro_values)
            )
        else:
            raise ValueError(f"Unsupported model input feature name: {feature}")

    matrix = pd.DataFrame(values, columns=required)
    if list(matrix.columns) != required:
        raise RuntimeError("Extracted feature columns do not match the model input schema.")
    if not np.isfinite(matrix.to_numpy(dtype=float)).all():
        raise ValueError("Extracted model feature matrix contains non-finite values.")
    return matrix


DEFAULT_MODEL_PATH = Path(__file__).resolve().parent.parent / "model" / "model.pkl"


def load_model(model_path: Path | None = None) -> object:
    """Load the packaged model artifact (default: static/model/model.pkl)."""
    path = model_path or DEFAULT_MODEL_PATH
    if not path.is_file():
        raise FileNotFoundError(
            f"Model file not found. Expected it at: {path}. "
            "Place model.pkl in the static/model folder."
        )
    with path.open("rb") as model_file:
        return pickle.load(model_file)
