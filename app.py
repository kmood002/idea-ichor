import streamlit as st
import pandas as pd

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)
st.subheader("Differential Expression Atlas (IDEA)")

@st.cache_data
def load_data():
    df = pd.read_excel("Murray_ProteinReport_26-118.xlsx", sheet_name="FullReport", header=1)
    df.columns = [str(col).strip() for col in df.columns]
    st.write("**Columns loaded:**", df.columns.tolist()[:15])   # Debug
    st.success(f"✅ Loaded {len(df):,} proteins")
    return df

df = load_data()

col1, col2, col3 = st.columns([3, 1.2, 1.2])
with col1:
    query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name", placeholder="Ca1, Alb, Gapdh")

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    st.write("**Searching for:**", terms)   # Debug

    mask = pd.Series(False, index=df.index)
    for col in df.columns:
        if any(k in col.lower() for k in ['gene', 'protein', 'name']):
            for term in terms:
                mask |= df[col].astype(str).str.upper().str.contains(term, na=False)
    
    res = df[mask].copy()
    st.success(f"Found {len(res)} rows before filtering")   # Debug

    # Rest of logic...
    fc_cols = [col for col in df.columns if 'DAY' in col and '/NAIVE' in col and not col.endswith('.1')]
    st.write("**FC columns:**", fc_cols)   # Debug

    if len(res) > 0:
        st.dataframe(res[['Genes', 'Protein Name']].head(10), use_container_width=True)
    else:
        st.error("No matches found — check debug info above")
