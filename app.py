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

query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name (comma-separated)", 
                     placeholder="Ca1, Alb, Gapdh, Col1a1, Actg1")

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    
    mask = pd.Series(False, index=df.index)
    for term in terms:
        mask |= df.iloc[:, 1].astype(str).str.upper().str.contains(term, na=False)   # Gene
        mask |= df.iloc[:, 2].astype(str).str.upper().str.contains(term, na=False)   # Protein Name
    
    res = df[mask].copy()
    
    # Identify Log2FC columns (they contain the fold changes)
    fc_cols = [col for col in df.columns if isinstance(col, str) and '/NAIVE' in col and not col.endswith('.1')]
    p_cols = [col for col in df.columns if isinstance(col, str) and col.endswith('.1') and 'NAIVE' in col]
    
    for c in fc_cols + p_cols:
        if c in res.columns:
            res[c] = pd.to_numeric(res[c], errors='coerce')
    
    if len(res) > 0 and len(fc_cols) >= 3:
        res['Max |log2FC|'] = res[fc_cols].abs().max(axis=1)
        res['Max Time'] = res[fc_cols].abs().idxmax(axis=1)
        
        def get_direction(row):
            mt = row['Max Time']
            val = row[mt] if mt in row else None
            return '↑ Up' if pd.notna(val) and val > 0 else '↓ Down'
        
        res['Direction'] = res.apply(get_direction, axis=1)
        
        # Clean display
        display = res.iloc[:, [1, 2]].copy()
        display.columns = ['Gene', 'Protein Name']
        display['Max |log2FC|'] = res['Max |log2FC|'].round(2)
        display['Max Time'] = res['Max Time']
        display['Direction'] = res['Direction']
        
        for i, col in enumerate(fc_cols):
            display[col.replace('/NAIVE', ' vs NAIVE (log2FC)')] = res[col].round(2)
        
        st.success(f"✅ Found {len(display)} matching proteins")
        st.dataframe(display, use_container_width=True, hide_index=True)
        st.download_button("📥 Download Results", res.to_csv(index=False), "dryeye_results.csv")
    else:
        st.warning("No matching genes found.")

st.caption("Basic search + fold-change display active. Ready for thresholds, multiple models, and further polish.")
