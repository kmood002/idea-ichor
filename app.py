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
    st.write("**Available columns:**", df.columns.tolist()[:30])  # Debug
    return df

df = load_data()

col1, col2, col3 = st.columns([3, 1.2, 1.2])
with col1:
    query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name", placeholder="Ca1, Alb, Gapdh")

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    st.write("**Searching terms:**", terms)
    
    mask = pd.Series(False, index=df.index)
    for col in df.columns:
        col_lower = str(col).lower()
        if any(k in col_lower for k in ['gene', 'protein', 'name']):
            for term in terms:
                mask |= df[col].astype(str).str.upper().str.contains(term, na=False)
    
    res = df[mask].copy()
    st.success(f"**Raw matches found:** {len(res)}")
    
    if len(res) > 0:
        st.dataframe(res.iloc[:, :6].head(5), use_container_width=True)  # Show raw data for debugging

    # ... (rest of filtering code can be added later)
