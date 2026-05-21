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

query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name (comma-separated)", placeholder="Ca13, Alb, Gapdh, Col1a1, Actg1")

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    
    mask = pd.Series(False, index=df.index)
    for term in terms:
        mask |= df.iloc[:, 1].astype(str).str.upper().str.contains(term, na=False)   # Gene column
        mask |= df.iloc[:, 2].astype(str).str.upper().str.contains(term, na=False)   # Protein Name column
    
    res = df[mask].copy()
    
    if len(res) > 0:
        st.success(f"✅ Found {len(res)} matching proteins")
        display = res.iloc[:, [1, 2]].copy()
        display.columns = ['Gene', 'Protein Name']
        st.dataframe(display, use_container_width=True, hide_index=True)
        st.download_button("📥 Download Matching Rows", res.to_csv(index=False), "matches.csv")
    else:
        st.warning("No matching genes found. Try 'Ca13', 'Alb', or 'Gapdh'")

st.caption("Basic search active. Thresholds and time-point highlighting will be re-added once search is confirmed working.")
