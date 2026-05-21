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

query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name", placeholder="Ca1, Alb, Gapdh, Col1a1, Actg1")

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    
    mask = pd.Series(False, index=df.index)
    for term in terms:
        mask |= df.iloc[:, 2].astype(str).str.upper().str.contains(term, na=False)   # Gene (index 2)
        mask |= df.iloc[:, 3].astype(str).str.upper().str.contains(term, na=False)   # Protein Name (index 3)
    
    res = df[mask].copy()
    
    if len(res) > 0:
        st.success(f"✅ Found {len(res)} matching proteins")
        
        # Clean display
        display = res.iloc[:, [2, 3]].copy()
        display.columns = ['Gene', 'Protein Name']
        
        # Add Log2FC columns (they start around column 23-25)
        fc_cols = df.columns[23:26].tolist()   # DAY2/NAIVE, DAY7/NAIVE, DAY14/NAIVE
        for i, col in enumerate(fc_cols):
            display[col] = pd.to_numeric(res[col], errors='coerce').round(2)
        
        # Max differential
        display['Max |log2FC|'] = display[fc_cols].abs().max(axis=1).round(2)
        display['Max Time'] = display[fc_cols].abs().idxmax(axis=1)
        
        st.dataframe(display, use_container_width=True, hide_index=True)
        st.download_button("📥 Download Results", res.to_csv(index=False), "dryeye_results.csv")
    else:
        st.warning("No matching genes found.")

st.caption("Test: Ca1, Alb, Gapdh, Col1a1")
