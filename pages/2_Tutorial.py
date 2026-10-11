from pathlib import Path
import streamlit as st

APP_ROOT = Path(__file__).resolve().parents[1]
STATIC_DIR = APP_ROOT / "static"
TUTORIAL_IMAGES = STATIC_DIR / "images" / "tutorial"

st.set_page_config(page_title="AMP-PepPI Tutorial", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image(str(STATIC_DIR / "images" / "icarlogo.png"), width=150)

with col2:
    st.markdown("<h2 style='text-align:center;'> AMP-PepPI: Antimicrobial Peptide-Protien Interaction tool</h2>", unsafe_allow_html=True)
st.markdown("---")

st.header("Tutorial on AMP-PepPI")

st.write(
    "This tutorial provides a step-by-step guide on predicting peptide-protein interactions (PPI) "
    "between antimicrobial peptide (AMP) and their target proteins using AMP-PepPI."
)
st.write(
    "Upon visiting the AMP-PepPI homepage, users can initiate the prediction process by clicking "
    "the “Prediction” button on the left side, as shown in Figure 1."
)
#st.image(caption="Figure 1: Homepage of AMP-PepPI.", image=TUTORIAL_IMAGES/"fig1_homepage.png",width=500)

st.subheader("Providing Input Sequences")
st.write("Users can input AMP and protein sequences in two ways:")
st.markdown(
    "- Paste the sequences directly into the respective text boxes.\n"
    "- Upload sequences in FASTA format."
)
st.write("Additionally, users can:")
st.markdown(
    "- Provide pre-defined AMP-protein pairs for interaction prediction.\n"
    "- Input a variable number of AMP and target proteins and generate all possible pairwise "
    "combinations for analysis (Figure 2)."
)

st.image(caption="a) Option to paste the sequence(s) in respective tab.", image=TUTORIAL_IMAGES/"fig2a_paste_input.png",width=500)
st.image(caption="b) Option to upload the sequence(s) in fasta format.", image=TUTORIAL_IMAGES/"fig2b_upload_all_combinations.png",width=500)
st.image(caption="c) Option to upload the sequence(s) in fasta format.", image=TUTORIAL_IMAGES/"fig2c_upload_paired.png",width=500)

st.caption("Figure 2: View of the prediction tab.")

st.subheader("Viewing and Interpreting Results")
st.write("Predictions are displayed in a tabular format and can be downloaded for further analysis (Figure 3).")
st.markdown("**Prediction outcomes:**")
st.markdown(
    "- 1 → Positive interaction predicted.\n"
    "- 0 → No interaction predicted."
)
st.write("Each prediction includes an interaction probability score.")

st.image(caption="Result table", image=TUTORIAL_IMAGES/"fig3a_results_table.png",width=500)
st.image(caption="Result charts", image=TUTORIAL_IMAGES/"fig3b_results_charts.png",width=500)

st.caption("Figure 3: Result view.")

st.write(
    "In the given example, the result shows that the probability of interaction between the first "
    "and third given protein pairs is 100%. ‘1’ in the ‘Prediction’ column means "
    "that the model predicted positive interaction. Whereas, second pair as well as the bottom five "
    "pairs have a higher probability of not interacting and hence the final prediction was in favor "
    "of negative interaction or no interaction represented by ‘0’."
)

st.subheader("Downloading and Analysing Results")
st.write(
    "Users can download results for offline analysis. The downloaded results can be viewed in Excel "
    "or any spreadsheet tool (Figure 4)."
)
st.image(caption="Figure 4: Snapshot of the result", image=TUTORIAL_IMAGES/"fig4_results_in_excel.png",width=500)

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)
