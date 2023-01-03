
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
import seaborn as sns 


st.set_page_config(layout="wide")


version = 0.1


    
#Set up main page of application / Header info / data collection / file selection / remove files / Reset

#Main landing page greating / info


st.title('Araya Deconvolution Setup Analysis Tool ' +str(version))

st.write("Upload parsed csv files for comparison - use ePCR Viewer to download 'all data csv' for each Araya before starting here - do this for each read set for your 'test Araya'")
st.write("Link to ePCR viewer [Here](https://share.streamlit.io/j0n0curry/new-calibration-test/main/ePCR_new_cal_OX_viewv1.py)")
st.write('Uploaded files - Reference is the Araya being used to compare to while test is the Araya being adjusted - no warranties supplied or implied for this application')
st.write('click on images to increase size')
####start loading data and set up for persistent dataframe via cache ans session state - TODO

@st.cache(allow_output_mutation=True)
def load_data(data):
    df = pd.read_csv(data)
    
    #df = normalise_values(df)
    return(df)

#####start of stateful - found streamlit blog - TODO - assign dropdown values and changes in data frame to select variables
##### dataframe will persistent through call session state and assigning the cvs to this - TODO - parser with caching for files#
##### after upload - likely pd.concat - reassign comp to new variable - cache here or allow mutation of orginal - horrible and messy 
##### rerun of small merge function to preseve and assign new df and then allow selection of columns or arrays etc....
def persist_dataframe():
    # drop column from dataframe
    delete_col = st.session_state["delete_col"]
    if delete_col in st.session_state["updated_df1"]:
        st.session_state["updated_df1"] = st.session_state["updated_df1"].drop(
            columns=[delete_col]
        )
    else:
        st.sidebar.warning("Column previously deleted. Select another column.")
    with col2:
        st.write("Updated dataframe")
        st.dataframe(st.session_state["updated_df1"])
        st.write(st.session_state["updated_df1"].columns.tolist())
  
######normalise values for z scores - absolute deviations from the mean - takes in all processes but could be used to call after session
######state is called to calculate abs deviations for a called group using result as a filter
def normalise_values(df): 
    df['norm_zscore'] = (df.ROX_RFU - df.ROX_RFU.mean())/df.ROX_RFU.std(ddof=0)
    df['cfo_zscore'] = (df.VIC_RFU - df.VIC_RFU.mean())/df.VIC_RFU.std(ddof=0)
    df['fam_zscore'] = (df.FAM_RFU - df.FAM_RFU.mean())/df.FAM_RFU.std(ddof=0)
    df['nFAM_zscore'] = (df.norm_N_Cov - df.norm_N_Cov.mean())/df.norm_N_Cov.std(ddof=0)
    df['nVIC_zscore'] = (df.norm_RNaseP - df.norm_RNaseP.mean())/df.norm_RNaseP.std(ddof=0)
    return(df)


######function to create percentage chage but avoids 0 0 div problems 
def pct_change(first, second):
        diff = second - first
        change = 0
        try:
            if diff > 0:
                change = (diff / first) * 100
            elif diff < 0:
                diff = first - second
                change = -((diff / first) * 100)
        except ZeroDivisionError:
            return float('inf')
        return change



uploaded_file1 = st.sidebar.file_uploader("Uploaded Reference Araya", type=['csv'], accept_multiple_files=False, key = 'key')


# style
th_props = [
  ('font-size', '30px'),
  ('text-align', 'center'),
  ('font-weight', 'bold'),
  ('color', '#6d6d6d'),
  ('background-color', '#f7ffff')
  ]
                               
td_props = [
  ('font-size', '12px')
  ]
                                 
styles = [
  dict(selector="th", props=th_props),
  dict(selector="td", props=td_props)
  ]





#load dataset and assign to dataframe - this is a session state frame - is static
if uploaded_file1 is not None:
    df1 = load_data(uploaded_file1)
else:
    st.warning('Please upload data for analysis')
    st.stop() 

# callback to session_state
# initialize session state variable
if "updated_df1" not in st.session_state:
    st.session_state.updated_df1 = df1
    


def apply_transformation(basis, T):
    """Return the transformed basis after applying transformation T."""
    return (T @ basis.T).T



s_f11 = 'F1/1'
s_f12 = 'F1/2'
s_f13 = 'F1/3'
s_f21 = 'F2/1'
s_f22 = 'F2/2'
s_f23 = 'F2/3'
s_f31 = 'F3/1'
s_f32 = 'F3/2'
s_f33 = 'F3/3'

a = st.sidebar.slider(s_f11,0.000, 2.000, 1.000)
b = st.sidebar.slider(s_f12, -1.000, 1.000, 0.000)
c = st.sidebar.slider(s_f13, -1.000, 1.000, 0.000)
d = st.sidebar.slider(s_f21, -1.000, 1.000, 0.000)
e = st.sidebar.slider(s_f22, 0.000, 2.000, 1.000)
f = st.sidebar.slider(s_f23, -1.000, 1.000, 0.000)
g = st.sidebar.slider(s_f31, -1.000, 1.000, 0.000)
h = st.sidebar.slider(s_f32, -1.000, 1.000, 0.000)
i = st.sidebar.slider(s_f33, 0.000, 2.000, 1.000)

st.write(a,b,c,d,e,f,g,h,i)
    
new_m = np.array([[a,b,c],[d,e,f],[g,h,i]])


adj_out = pd.DataFrame(new_m, columns = [['FAM_IN', 'VIC_IN', 'ROX_IN']],index = [['FAM_out', 'VIC_out', 'ROX_out']])

    # table
df2 = adj_out.style.set_properties(**{'text-align': 'left'}).set_table_styles(styles)
st.table(df2)




basis = df1[['FAM_RFU', 'VIC_RFU', 'ROX_RFU']].to_numpy()
S = new_m
tbasis = apply_transformation(basis, S)

df1[['new_FAM', 'new_VIC','new_ROX']] = tbasis

df1['comp_nFAM'] = df1.new_FAM / df1.new_ROX

df1['comp_nVIC'] = df1.new_VIC / df1.new_ROX



#vals, vecs = np.linalg.eig(S)



#fig = plt.gcf()
#fig.set_size_inches(20,12)

def plot_raw_over(df,a,b,c,d,L,W):
    fig1, (ax_nFAM, ax_fam_rox) = plt.subplots(nrows=1,ncols=2,
    figsize=(L,W))
    
    sns.scatterplot(data = df1, x = a, y = b, s = 30, color = 'red', edgecolor = 'black', ax=ax_nFAM)
    

    sns.scatterplot(data = df1, x = c, y = d, s = 30, hue='Result', edgecolor = 'black', ax=ax_nFAM)
    
    sns.scatterplot(data = df1,x= 'new_ROX', y = 'comp_nFAM', s = 30, color = 'red', edgecolor = 'black', ax=ax_fam_rox)
    
    
    ax_nFAM.grid(True)
    ax_fam_rox.grid(True)
    st.pyplot(fig1)

plot_raw_over(df1,'comp_nVIC','comp_nFAM','norm_RNAseP','norm_N_Cov', 20,10)





#st.pyplot(sns.scatterplot(data = df1,x= 'comp_nVIC', y = 'comp_nFAM', s = 30, color = 'red', edgecolor = 'black'))
#sns.scatterplot(data = df1,x= 'norm_RNAseP', y = 'norm_N_Cov', s = 30, color = 'blue', edgecolor = 'black')

#st.pyplot(plt)


#sns.scatterplot(data = df1,x= 'new_ROX', y = 'comp_nFAM', s = 5, color = 'red', edgecolor = 'black')
#plt.grid()  #just add this
#st.pyplot(plt)
   
 
