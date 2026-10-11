import streamlit as st

st.set_page_config(page_title="AMP-PepPI Tutorial", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image("static/images/icarlogo.png", width=150)

with col2:
    st.markdown("<h2 style='text-align:center;'> AMP-PepPI: Antimicrobial Peptide-Protein Interaction tool</h2>", unsafe_allow_html=True)

st.markdown("---")
st.text("")

col1_1, col2_1 = st.columns([1, 1])

with col1_1:
    with st.container(border=True, height=160):
        st.markdown('''<div style='text-align:center'><b>Deeksha P.M.</b></br>  
                    Discipline of Bioinformatics, The Graduate School  </br>
                    ICAR-Indian Agricultural Research Institute  </br>
                    Pusa, New Delhi-110012  </br>
                    Mail Id: 27deekshapm@gmail.com</div>''', unsafe_allow_html=True)

with col2_1:
     with st.container(border=True, height=160):
        st.markdown('''<div style='text-align:center'><b>Dr. Shashi Bhusan Lal</b></br>
                    Principal Scientist, ICAT-Unit </br>  
                    Indian Council of Agricultural Research </br>
                    Pusa, New Delhi-110012  </br>
                    Mail Id: sb.lal@icar.org.in</div>''', unsafe_allow_html=True)


st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)
