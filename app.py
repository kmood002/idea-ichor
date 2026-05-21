import streamlit as st
import pandas as pd

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)
st.subheader("Differential Expression Atlas (IDEA)")
st.caption("**Model:** Scopolamine + Desiccating Stress Dry Eye | **Tissue:** Cornea | C57BL/6 Mice")

@st.cache_data
def load_data():
    df = pd.read_excel("Murray_ProteinReport_26-118.xlsx", sheet_name="FullReport", header=2)
    df.columns = [str(col).strip() for col in df.columns]
    st.success(f"✅ Loaded {len(df):,} proteins")
    return df

df = load_data()

st.write("**First few rows for debugging:**")
st.dataframe(df.iloc[:, :8].head(3), use_container_width=True)

col1, col2, col3 = st.columns([3, 1.2, 1.2])
with col1:
    query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name", placeholder="Ca1, Alb, Gapdh")

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    st.write("**Searching for:**", terms)
    
    mask = pd.Series(False, index=df.index)
    for term in terms:
        mask |= df.iloc[:, 1].astype(str).str.upper().str.contains(term, na=False)  # Gene column
        mask |= df.iloc[:, 2].astype(str).str.upper().str.contains(term, na=False)  # Protein Name column
    
    st.write(f"**Raw matches found:** {mask.sum()}")
    
    res = df[mask].copy()
    
    if len(res) > 0:
        st.dataframe(res.iloc[:, :8], use_container_width=True)
    else:
        st.error("No matches — please tell me what the debug info above shows.")
