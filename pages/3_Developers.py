import streamlit as st

st.set_page_config( page_title="AMP-PepPI", initial_sidebar_state="expanded", layout="wide")
col1, col2 = st.columns([1.5, 20])

with col1:
    st.image("static/images/icarlogo.png", width=150)

with col2:
    st.markdown("<h1 style='text-align:center;'> AMP-PepPI: A Machine-Learning based Tool for Predicting Interactions between Antimicrobial Peptides and Target Proteins</h1>", unsafe_allow_html=True)


st.markdown("---")
st.text("")

col1_1, col2_1 =st.columns([1, 1])

with col1_1:
    with st.container(border=True):
        colInCon_1, colInCon_2 = st.columns([1, 3])
        with colInCon_1:
            pass#st.image("static/images/deeksha.jpg")
        with colInCon_2:
            st.markdown('''**Deeksha P.M.**  
                        Ph.D. Scholar (Bioinformatics)  
                        ICAR-Indian Agricultural Statistics Research Institute  
                        Pusa, New Delhi-110012, India.  
                        Contact mail: 27deekshapm@gmail.com''')
    
    with st.container(border=True):
        colInCon_1, colInCon_2 = st.columns([1, 3])
        with colInCon_1:
            pass#st.image("static/images/sneha.jpg")
        with colInCon_2:
            st.markdown('''**Dr. Sneha Murmu**  
                        Scientist  
                        ICAR-Indian Agricultural Statistics Research Institute  
                        Pusa, New Delhi-110012, India.  
                        Contact mail: murmu.sneha07@gmail.com''')
            
    with st.container(border=True):
        colInCon_1, colInCon_2 = st.columns([1, 3])
        with colInCon_1:
            pass#st.image("static/images/dwijesh.jpeg")
        with colInCon_2:
            st.markdown('''**Dr. Dwijesh Chandra Mishra**  
                        Senior Scientist  
                        ICAR-Indian Agricultural Statistics Research Institute  
                        Pusa, New Delhi-110012, India.     
                        Contact mail: dwij.mishra@gmail.com''')
            

with col2_1:
    with st.container(border=True):
        colInCon_1, colInCon_2 = st.columns([1, 3])
        with colInCon_1:
            pass#st.image("static/images/sblal.jpg")
        with colInCon_2:
            st.markdown('''**Dr. Shashi Bhushan Lal**  
                        Principal Scientist  
                        ICT Unit, Indian Council of Agricultural Research  
                        Pusa, New Delhi-110012, India.  
                        Contact mail: sb.lal@icar.org.in''')
            
    with st.container(border=True):
        colInCon_1, colInCon_2 = st.columns([1, 3])
        with colInCon_1:
            pass#st.image("static/images/himanshu.png")
        with colInCon_2:
            st.markdown('''**Dr. Himanshushekhar Chaurasia**  
                        Scientist (Computer Application & IT)  
                        ICAR - Central Institute for Research on Cotton Technology,  
                        Mumbai, Maharashtra-400019, India.  
                        Contact mail: h.chaurasia@icar.org.in''')
            
    with st.container(border=True):
        colInCon_1, colInCon_2 = st.columns([1, 3])
        with colInCon_1:
            pass#st.image("static/images/joshitha.jpg")
        with colInCon_2:
            st.markdown('''**Dr. Joshitha Vijayan**  
                        Scientist  
                        ICAR-Nantional Institute for Plant Biotechnology  
                        New Delhi-110012, India.  
                        Contact mail: jovijayan@gmail.com''')
            

st.text("")
st.markdown("<div style='background-color:#32CD32; text-align:center'><p style='color:white'>Copyright © 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi-110012. All rights reserved.</p></div>", unsafe_allow_html=True)
