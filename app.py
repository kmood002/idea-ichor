import streamlit as st
import pandas as pd

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)
st.subheader("Differential Expression Atlas (IDEA)")
st.caption("**Model:** Scopolamine + Desiccating Stress Dry Eye | **Tissue:** Cornea | C57BL/6 Mice")

@st.cache_data
def load_data():
    df = pd.read_excel("Murray_ProteinReport_26-118.xlsx", sheet_name="FullReport", header=1)
    df.columns = [str(col).strip() for col in df.columns]
    st.success(f"✅ Loaded {len(df):,} proteins")
    return df

df = load_data()

# Search & Filters
col1, col2, col3 = st.columns([3, 1.2, 1.2])
with col1:
    query = st.text_input("🔍 Enter Gene Symbol(s) or Protein Name (comma-separated)", placeholder="Alb, Gapdh, Col1a1, Actg1")
with col2:
    fc_thresh = st.slider("|log2FC| ≥", 0.0, 5.0, 1.0, 0.1)
with col3:
    p_thresh = st.slider("p-value <", 0.0001, 0.1, 0.05, 0.001)

if query:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    
    # Robust column search
    gene_col = next((col for col in df.columns if 'gene' in col.lower()), 'Genes')
    name_col = next((col for col in df.columns if 'protein name' in col.lower() or 'name' in col.lower()), 'Protein Name')
    
    mask = pd.Series([False] * len(df))
    if gene_col in df.columns:
        mask |= df[gene_col].astype(str).str.upper().str.contains('|'.join(terms), na=False)
    if name_col in df.columns:
        mask |= df[name_col].astype(str).str.upper().str.contains('|'.join(terms), na=False)
    
    res = df[mask].copy()
    
    fc_cols = [col for col in df.columns if col.startswith('DAY') and '/NAIVE' in col and not col.endswith('.1')]
    p_cols = [col for col in df.columns if col.startswith('DAY') and '/NAIVE' in col and col.endswith('.1')]
    
    for c in fc_cols + p_cols:
        if c in res.columns:
            res[c] = pd.to_numeric(res[c], errors='coerce')
    
    res['Max |log2FC|'] = res[fc_cols].abs().max(axis=1)
    res['Max Time'] = res[fc_cols].abs().idxmax(axis=1)
    res['Direction'] = res.apply(lambda r: '↑ Up' if pd.notna(r[r['Max Time']]) and r[r['Max Time']] > 0 else '↓ Down', axis=1)
    
    filtered = res[(res['Max |log2FC|'] >= fc_thresh) & (res[p_cols].min(axis=1) < p_thresh)]
    
    if not filtered.empty:
        st.success(f"✅ Found {len(filtered)} matching proteins")
        display_cols = ['Genes', 'Protein Name', 'Max |log2FC|', 'Max Time', 'Direction'] + fc_cols
        st.dataframe(filtered[display_cols].head(50), use_container_width=True, hide_index=True)  # limit for performance
        st.download_button("📥 Download Full Results CSV", filtered.to_csv(index=False), "dryeye_results.csv")
    else:
        st.warning("No proteins meet the current thresholds.")

st.caption("Initial cornea dataset loaded. Ready for expansion.")
