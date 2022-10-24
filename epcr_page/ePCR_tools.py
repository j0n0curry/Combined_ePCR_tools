import streamlit as st

#to add a new page, simply add in the new appp in to the 'pages' folder and add in a function to call the new app below.
#Use the function to create a dropdown selectable for the app, place the list of apps and he markdown text to be show in the 
# pages_names_to_funcs dictionary. This will present the new app. 

def main_page():
    st.markdown("# ePCR tools : A selection of tools for calling, viewing an analysing data for ePCR and Arayas. Each page will have individual instruction to follow. If you need assistance email jonathan.curry@lgcgroup.com")
    st.sidebar.markdown("# ePCR tools 🎈")
    
def ePCR():
    st.markdown("# ePCR viewer ❄️")
    st.sidebar.markdown("# ePCR viewer ❄️")
    
def Araya():
    st.markdown("# Araya calibration 🎉")
    st.sidebar.markdown("# Araya verification 🎉")
    
def ePCR_RFL():
    st.markdown("# Araya calibration 🎉")
    st.sidebar.markdown("# Araya verification 🎉")
    
    
page_names_to_funcs = {
    "Main Page": main_page,
    "ePCR_viewer": ePCR,
    "Araya_comparison": Araya,
    "ePCR_Rosalind_Franklin": ePCR_RFL
}

selected_page = st.sidebar.selectbox("Select a page", page_names_to_funcs.keys())
page_names_to_funcs[selected_page]()