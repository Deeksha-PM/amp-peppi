import streamlit as st

st.set_page_config(page_title="AMP-PepPI Tutorial", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image("static/images/icarlogo.png", width=150)

with col2:
    st.markdown("<h2 style='text-align:center;'> AMP-PepPI: Antimicrobial Peptide-Protein Interaction tool</h2>", unsafe_allow_html=True)



st.markdown("---")
st.text("")

col1_1, col2_1 = st.columns([2, 2])

with col1_1:
    st.header("Background")
    st.text("Antimicrobial peptides (AMPs) are important components of innate immunity in diverse organisms. Understanding how AMPs interact with target proteins can provide insights into molecular recognition and support the identification of candidate peptide-protein interactions for further investigation.\nAMP-PepPI is a sequence-based machine learning tool developed to predict peptide-protein interactions (PepPIs) involving antimicrobial peptides and their target proteins. The predictor uses an Extreme Gradient Boosting (XGBoost) classifier trained on a combination of three complementary sequence-derived feature descriptors: Amino Acid Composition (AAC), Composition/Transition/Distribution (CTD), and physicochemical properties.\nAMP-PepPI provides a computational approach for prioritizing potential AMP–protein interactions using amino acid sequence information.The predictions are intended to support hypothesis generation and the prioritization of candidate interactions for subsequent experimental validation.", text_alignment="justify")
    if st.button("Start Prediction", type="primary"):
        st.switch_page("pages/1_Prediction.py")
with col2_1:
    st.image("static/images/Workflow.png")

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)
