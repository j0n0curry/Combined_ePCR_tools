import streamlit as st
import os
from io import BytesIO
from io import StringIO
from io import TextIOWrapper
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import math



def main_page():
    st.markdown("# ePCR tools")
    st.markdown("Providing tools for analysis and calling of data from Araya 1 / 2 - alpha state. The Araya comparison tool is intended to function with calibration tape and not to compare Arayas directly. Please down load parsed data - all_data.csv to perform separate analysis.")
    st.markdown("For question or improvement - please email jonathan.curry@lgcgroup.com with files or suggestions')
    st.sidebar.markdown("# Main page ")
    
def Araya1():
    st.markdown("# ePCR viewer - Araya 1 - 100 100 100")
    st.sidebar.markdown("# ePCR viewer - Araya 1")
    
def Araya_comp():
    st.markdown("# Araya calibration ")
    st.sidebar.markdown("# Araya verification ")
    
def Araya2():
    st.markdown("# ePCR viewer - Araya 2 - 100 100 100")
    st.sidebar.markdown("# ePCR viewer - Araya 2")
    
    
page_names_to_funcs = {
    "Main Page": main_page,
    "ePCR_viewer": Araya1,
    "Araya_comparison": Araya_comp,
    "ePCR_viewer" : Araya2
}

selected_page = st.sidebar.selectbox("Select a page", page_names_to_funcs.keys())
page_names_to_funcs[selected_page]()
