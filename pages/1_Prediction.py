import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

APP_ROOT = Path(__file__).resolve().parents[1]
STATIC_DIR = APP_ROOT / "static"
SCRIPT_DIR = STATIC_DIR / "script"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from feature_extraction import (  # noqa: E402
    AMINO_ACIDS,
    AMP_MAX_LEN,
    AMP_MIN_LEN,
    PROTEIN_MAX_LEN,
    PROTEIN_MIN_LEN,
    load_model,
)
from predict import predict_frame  # noqa: E402
from validation import (  # noqa: E402
    MAX_BULK_PAIRS,
    count_fasta_records,
    length_status,
    parse_fasta,
    parse_single_input,
    validate_records,
)

EXAMPLE_AMP = STATIC_DIR / "input" / "example_AMPs.fasta"
EXAMPLE_PROTEIN = STATIC_DIR / "input" / "example_proteins.fasta"
FASTA_TYPES = ["fasta", "fa", "faa", "fas", "txt"]
CHUNK_SIZE = 500
MAX_SHOWN_MESSAGES = 10
STATE_KEYS = ("results", "errors", "invalid", "warnings")

st.set_page_config(page_title="AMP-PepPI", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image(str(STATIC_DIR / "images" / "icarlogo.png"), width=150)

with col2:
    st.markdown("<h2 style='text-align:center;'> AMP-PepPI: Antimicrobial Peptide-Protien Interaction tool</h2>", unsafe_allow_html=True)
st.markdown("---")
st.text("")


@st.cache_resource
def get_model():
    return load_model()


def read_text(path):
    return Path(path).read_text()


def read_upload(uploaded_file):
    return uploaded_file.getvalue().decode("utf-8", errors="replace")


def safe_record_count(text):
    try:
        return count_fasta_records(text)
    except Exception:
        return 0


def reset_outputs():
    for key in STATE_KEYS:
        st.session_state[key] = None


def run_prediction(pairs):
    """pairs: list of (amp_id, amp_seq, protein_id, protein_seq). Returns a results frame or None."""
    progress = st.progress(0.0, text="Calculating features and predicting...")
    try:
        model = get_model()
        frame = pd.DataFrame(
            {
                "amp": [p[0] for p in pairs],
                "protein": [p[2] for p in pairs],
                "amp_sequence": [p[1] for p in pairs],
                "protein_sequence": [p[3] for p in pairs],
            }
        )

        def update(done, total):
            progress.progress(done / total, text=f"Processed {done} of {total} pairs")

        predictions = predict_frame(frame, model, chunk_size=CHUNK_SIZE, progress_callback=update)
        results = pd.concat([frame, predictions], axis=1)
        results.insert(2, "Prediction", ["Positive" if y == 1 else "Negative" for y in results["Y_pred"]])
        return results
    except FileNotFoundError as e:
        st.session_state["errors"] = [str(e)]
    except Exception as e:
        st.session_state["errors"] = [f"Prediction error: {e}"]
    finally:
        progress.empty()
    return None


def process(amp_records, prot_records, mode):
    """Validate the records, build the pairs ("all" = all-vs-all, "paired" = 1-to-1) and predict."""
    reset_outputs()
    amps = validate_records(amp_records, "AMP")
    prots = validate_records(prot_records, "Protein")
    st.session_state["warnings"] = amps.warnings + prots.warnings
    errors = amps.errors + prots.errors

    if errors:
        st.session_state["errors"] = errors
        st.session_state["invalid"] = amps.invalid + prots.invalid
        return

    if mode == "all":
        n_pairs = len(amps.valid) * len(prots.valid)
        if n_pairs > MAX_BULK_PAIRS:
            st.session_state["errors"] = [
                f"{len(amps.valid)} AMPs x {len(prots.valid)} proteins = {n_pairs} pairs, which is "
                f"more than the maximum of {MAX_BULK_PAIRS} pairs per run. Please split the input."
            ]
            return
        pairs = [(a["id"], a["seq"], p["id"], p["seq"]) for a in amps.valid for p in prots.valid]
    else:
        n = min(len(amps.valid), len(prots.valid))
        pairs = [
            (amps.valid[i]["id"], amps.valid[i]["seq"], prots.valid[i]["id"], prots.valid[i]["seq"])
            for i in range(n)
        ]
        if len(amps.valid) != len(prots.valid):
            st.session_state["warnings"].append(
                f"File lengths differ ({len(amps.valid)} AMPs, {len(prots.valid)} proteins). "
                f"Only the first {n} pairs were processed."
            )

    st.session_state["results"] = run_prediction(pairs)


def process_fasta(amp_text, prot_text, mode):
    errors = []
    try:
        amp_records = parse_fasta(amp_text)
    except ValueError as e:
        errors.append(f"AMP file: {e}")
    try:
        prot_records = parse_fasta(prot_text)
    except ValueError as e:
        errors.append(f"Protein file: {e}")
    if errors:
        reset_outputs()
        st.session_state["errors"] = errors
        return
    process(amp_records, prot_records, mode)


def render_messages(messages, show):
    """Show the first few messages individually and fold the rest into an expander."""
    for msg in messages[:MAX_SHOWN_MESSAGES]:
        show(msg)
    rest = messages[MAX_SHOWN_MESSAGES:]
    if rest:
        with st.expander(f"Show {len(rest)} more"):
            st.markdown("\n".join(f"- {m}" for m in rest))


st.header("Interaction Prediction")

st.info(
    f"Input requirements: AMP length {AMP_MIN_LEN}-{AMP_MAX_LEN} residues, protein length "
    f"{PROTEIN_MIN_LEN}-{PROTEIN_MAX_LEN} residues. Only the 20 standard amino acids "
    f"({AMINO_ACIDS}) are accepted. Download the example files to see how the data should be formatted."
)
c_dl1, c_dl2 = st.columns(2)
with c_dl1:
    if EXAMPLE_AMP.exists():
        st.download_button("Download Example AMPs (FASTA)", EXAMPLE_AMP.read_bytes(),
                           "example_amps.fasta", on_click="ignore")
with c_dl2:
    if EXAMPLE_PROTEIN.exists():
        st.download_button("Download Example Proteins (FASTA)", EXAMPLE_PROTEIN.read_bytes(),
                           "example_proteins.fasta", on_click="ignore")

st.markdown("---")

tab1, tab2, tab3 = st.tabs(["Single Pair", "All Combinations", "Paired Upload"])

with tab1:
    st.subheader("Paste the Sequences")
    st.caption("Paste a plain sequence or a FASTA record (only the first record is used). "
               "The length counter updates when you click outside the box or press Ctrl+Enter.")
    c1, c2 = st.columns(2)
    with c1:
        s_amp = st.text_area("AMP Sequence", placeholder="KLSPSL...")
        st.caption(length_status(s_amp, "AMP"))
    with c2:
        s_prot = st.text_area("Protein Sequence", placeholder="MRPPQC...")
        st.caption(length_status(s_prot, "Protein"))

    if st.button("Predict Single Interaction", type="primary"):
        try:
            amp_id, amp_raw, amp_warn = parse_single_input(s_amp, "Manual_AMP")
            prot_id, prot_raw, prot_warn = parse_single_input(s_prot, "Manual_Protein")
        except ValueError as e:
            reset_outputs()
            st.session_state["errors"] = [str(e)]
        else:
            process([(amp_id, amp_raw)], [(prot_id, prot_raw)], "paired")
            if st.session_state["warnings"] is not None:
                st.session_state["warnings"] = amp_warn + prot_warn + st.session_state["warnings"]

with tab2:
    st.subheader("All Combinations (All-vs-All)")
    u1, u2 = st.columns(2)
    with u1:
        file_amp = st.file_uploader("Upload amp.fasta", type=FASTA_TYPES, key="bulk_u_amp")
    with u2:
        file_prot = st.file_uploader("Upload protein.fasta", type=FASTA_TYPES, key="bulk_u_prot")

    use_example_bulk = st.checkbox("Use example sequences instead")

    amp_text = prot_text = None
    if use_example_bulk and EXAMPLE_AMP.exists() and EXAMPLE_PROTEIN.exists():
        amp_text, prot_text = read_text(EXAMPLE_AMP), read_text(EXAMPLE_PROTEIN)
    elif file_amp and file_prot:
        amp_text, prot_text = read_upload(file_amp), read_upload(file_prot)

    if amp_text is not None:
        n_a, n_p = safe_record_count(amp_text), safe_record_count(prot_text)
        st.caption(f"{n_a} AMPs x {n_p} proteins = {n_a * n_p} pairs (maximum {MAX_BULK_PAIRS} per run).")

    if st.button("Run Bulk Analysis", type="primary"):
        if amp_text is not None:
            process_fasta(amp_text, prot_text, "all")
        else:
            st.warning("Please upload both files or tick 'Use example sequences instead'.")

with tab3:
    st.subheader("Paired Prediction (1-to-1)")
    st.caption("Upload separate files to predict specific pairs (AMP 1 with Protein 1, etc.) instead of all combinations.")

    u3_col1, u3_col2 = st.columns(2)
    with u3_col1:
        file_amp_p = st.file_uploader("Upload amp.fasta", type=FASTA_TYPES, key="paired_u_amp")
    with u3_col2:
        file_prot_p = st.file_uploader("Upload protein.fasta", type=FASTA_TYPES, key="paired_u_prot")

    run_example_paired = st.button("Run Example Paired Prediction")

    if st.button("Predict Pairs", type="primary") or run_example_paired:
        if run_example_paired:
            if EXAMPLE_AMP.exists() and EXAMPLE_PROTEIN.exists():
                process_fasta(read_text(EXAMPLE_AMP), read_text(EXAMPLE_PROTEIN), "paired")
            else:
                st.warning("Example files were not found in static/data.")
        elif file_amp_p and file_prot_p:
            process_fasta(read_upload(file_amp_p), read_upload(file_prot_p), "paired")
        else:
            st.warning("Please upload both files or click 'Run Example Paired Prediction'.")

# Outputs (kept in session state so they survive reruns caused by the filter widgets)
errors = st.session_state.get("errors")
invalid = st.session_state.get("invalid")
warnings = st.session_state.get("warnings")
results = st.session_state.get("results")

if errors or warnings:
    st.markdown("---")
if errors:
    st.subheader("Input errors")
    st.caption("No prediction was run. Fix the problems below and try again.")
    render_messages(errors, st.error)
    if invalid:
        st.download_button(
            "Download invalid records (CSV)",
            pd.DataFrame(invalid).to_csv(index=False).encode("utf-8"),
            "invalid_records.csv", "text/csv", on_click="ignore",
        )
if warnings:
    render_messages(warnings, st.warning)

if results is not None:
    st.markdown("---")
    st.subheader("Results & Analytics")

    PROB = "Interaction probability"
    results = results.assign(**{PROB: results["probability_class_1"]})

    if len(results) == 1:
        st.metric("Predicted interaction probability", f"{results[PROB].iloc[0]:.4f}")
        st.success(f"Result: {results['Prediction'].iloc[0]}")

    f1, f2, f3, f4 = st.columns(4)
    with f1:
        min_prob = st.slider("Minimum interaction probability", 0.0, 1.0, 0.0, 0.01)
    with f2:
        label_filter = st.selectbox("Prediction", ["All", "Positive", "Negative"])
    with f3:
        top_n = st.number_input("Top N proteins per AMP (0 = all)", min_value=0, value=0, step=1)
    with f4:
        sort_by = st.selectbox(
            "Sort by", ["Original order", "Probability (high to low)", "Probability (low to high)"]
        )

    view = results[results[PROB] >= min_prob]
    if label_filter != "All":
        view = view[view["Prediction"] == label_filter]
    if sort_by != "Original order" or top_n > 0:
        view = view.sort_values(PROB, ascending=(sort_by == "Probability (low to high)"))
    if top_n > 0:
        view = view.groupby("amp", sort=False).head(int(top_n))
    if sort_by == "Original order":
        view = view.sort_index()

    st.caption(f"Showing {len(view)} of {len(results)} predictions.")

    if view.empty:
        st.info("No predictions match the current filters.")
    else:
        c_g1, c_g2 = st.columns([2, 1])
        with c_g1:
            fig = px.histogram(
                view, x=PROB, color="Prediction",
                title="Interaction Probability Distribution",
                color_discrete_map={"Positive": "#00CC96", "Negative": "#EF553B"},
                marginal="box",
            )
            st.plotly_chart(fig, width="stretch")

        with c_g2:
            pie_fig = px.pie(
                view, names="Prediction", title="Summary of Interactions",
                color="Prediction", color_discrete_map={"Positive": "#00CC96", "Negative": "#EF553B"},
            )
            st.plotly_chart(pie_fig, width="stretch")

        shown = view[["amp", "protein", "Prediction", PROB]].rename(
            columns={"amp": "AMP", "protein": "Protein"}
        )
        st.dataframe(shown.style.background_gradient(subset=[PROB], cmap="RdYlGn"), width="stretch")

    export_cols = [c for c in results.columns if c != PROB]
    d1, d2 = st.columns(2)
    with d1:
        st.download_button(
            "Download filtered results (CSV)",
            view[export_cols].to_csv(index=False).encode("utf-8"),
            "results_filtered.csv", "text/csv", on_click="ignore",
        )
    with d2:
        st.download_button(
            "Download all results (CSV)",
            results[export_cols].to_csv(index=False).encode("utf-8"),
            "results_all.csv", "text/csv", on_click="ignore",
        )

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)
