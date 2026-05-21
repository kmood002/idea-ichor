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

col1, col2, col3 = st.columns([3, 1.2, 1.2])
with col1:
    query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name", placeholder="Ca1, Alb, Gapdh, Col1a1, Actg1")
with col2:
    fc_thresh = st.slider("|log2FC| ≥", 0.0, 5.0, 1.0, 0.1)
with col3:
    p_thresh = st.slider("p-value <", 0.0001, 0.1, 0.05, 0.001)

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    
    mask = pd.Series(False, index=df.index)
    for term in terms:
        mask |= df.iloc[:, 1].astype(str).str.upper().str.contains(term, na=False)
        mask |= df.iloc[:, 2].astype(str).str.upper().str.contains(term, na=False)
    
    res = df[mask].copy()
    
    fc_cols = df.columns[23:26].tolist()   # DAY2/NAIVE, DAY7/NAIVE, DAY14/NAIVE
    p_cols = df.columns[26:29].tolist()
    
    for c in fc_cols + p_cols:
        if c in res.columns:
            res[c] = pd.to_numeric(res[c], errors='coerce')
    
    if len(res) > 0 and len(fc_cols) > 0:
        res['Max |log2FC|'] = res[fc_cols].abs().max(axis=1)
        res['Max Time'] = res[fc_cols].abs().idxmax(axis=1)
        
        def get_direction(row):
            mt = row['Max Time']
            val = row[mt] if mt in row else None
            return '↑ Up' if pd.notna(val) and val > 0 else '↓ Down'
        
        res['Direction'] = res.apply(get_direction, axis=1)
        
        filtered = res[(res['Max |log2FC|'] >= fc_thresh) & (res[p_cols].min(axis=1, skipna=True) < p_thresh)]
        
        if not filtered.empty:
            st.success(f"✅ Found {len(filtered)} matching proteins")
            display = filtered.iloc[:, [1, 2]].copy()  # Gene, Protein Name
            display.columns = ['Gene', 'Protein Name']
            display['Max |log2FC|'] = filtered['Max |log2FC|']
            display['Max Time'] = filtered['Max Time']
            display['Direction'] = filtered['Direction']
            display = pd.concat([display, filtered[fc_cols]], axis=1)
            st.dataframe(display, use_container_width=True, hide_index=True)
            st.download_button("📥 Download Results", filtered.to_csv(index=False), "dryeye_results.csv")
        else:
            st.warning("No proteins meet the thresholds.")
    else:
        st.info("No matching genes found.")

st.caption("Ready for additional models/tissues. Let me know what to improve next.")
