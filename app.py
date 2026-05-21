import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

# Logo
logo_path = "logo.png"
if Path(logo_path).exists():
    st.image(logo_path, width=320)
else:
    st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)

st.subheader("Ichor Differential Expression Atlas (IDEA)")
st.caption("**Pre-Clinical Model Explorer**")

# Load all data-*.xlsx files
data_files = list(Path(".").glob("data-*.xlsx"))

models = {}
for file in data_files:
    model_key = file.stem.replace("data-", "")
    xls = pd.ExcelFile(file)
    
    cover = pd.read_excel(xls, "cover", header=None)
    meta = dict(zip(cover.iloc[:, 0].astype(str).str.strip(), cover.iloc[:, 1].astype(str).str.strip()))
    
    model_name = meta.get("Model", model_key) if "Model" in meta else (cover.iloc[0, 1] if len(cover) > 0 else model_key)
    
    log2 = pd.read_excel(xls, "log2")
    pval = pd.read_excel(xls, "p")
    
    models[model_name] = {
        "meta": meta,
        "log2": log2,
        "pval": pval
    }

# Sidebar Filters
st.sidebar.header("Filters")
fc_thresh = st.sidebar.slider("|log2FC| ≥ (0 = show all)", 0.0, 5.0, 0.0, 0.1)
p_thresh = st.sidebar.slider("p-value < (0.1 = show all)", 0.0001, 0.1, 0.1, 0.001)

# Sidebar Column Toggles
st.sidebar.header("Display Columns")
show_model = st.sidebar.checkbox("Model", value=True)
show_gene = st.sidebar.checkbox("Gene", value=True)
show_protein = st.sidebar.checkbox("Protein Name", value=True)
show_species = st.sidebar.checkbox("Species", value=True)
show_strain = st.sidebar.checkbox("Strain", value=True)
show_gender = st.sidebar.checkbox("Gender", value=True)
show_tissue = st.sidebar.checkbox("Tissue", value=True)
show_pvalues = st.sidebar.checkbox("P-values", value=True)

# Search
query = st.text_input("🔍 Search by Protein Accession, Gene, or Protein Name", 
                     placeholder="Ca1, Alb, Gapdh, Col1a1")

if query and models:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    results = []
    
    for model_name, data in models.items():
        log2_df = data["log2"]
        p_df = data["pval"]
        meta = data["meta"]
        
        mask = (log2_df.iloc[:, 0].astype(str).str.upper().str.contains('|'.join(terms), na=False)) | \
               (log2_df.iloc[:, 1].astype(str).str.upper().str.contains('|'.join(terms), na=False)) | \
               (log2_df.iloc[:, 2].astype(str).str.upper().str.contains('|'.join(terms), na=False))
        
        hits = log2_df[mask].copy()
        if len(hits) > 0:
            for idx, row in hits.iterrows():
                d2_fc = row.iloc[3] if len(row) > 3 else None
                d7_fc = row.iloc[4] if len(row) > 4 else None
                d14_fc = row.iloc[5] if len(row) > 5 else None
                
                d2_p = p_df.iloc[idx, 3] if len(p_df.columns) > 3 else None
                d7_p = p_df.iloc[idx, 4] if len(p_df.columns) > 4 else None
                d14_p = p_df.iloc[idx, 5] if len(p_df.columns) > 5 else None
                
                include = True
                if fc_thresh > 0:
                    include &= any(abs(x or 0) >= fc_thresh for x in [d2_fc, d7_fc, d14_fc])
                if p_thresh < 0.1:
                    include &= any((x or 1) < p_thresh for x in [d2_p, d7_p, d14_p])
                
                if include:
                    results.append({
                        "Protein Accession": row.iloc[0],
                        "Gene": row.iloc[1],
                        "Protein Name": row.iloc[2],
                        "Model": model_name,
                        "Species": meta.get("Species", ""),
                        "Strain": meta.get("Strain", ""),
                        "Gender": meta.get("Gender", ""),
                        "Tissue": meta.get("Tissue", ""),
                        "Day 2 log2FC": round(d2_fc, 2) if d2_fc is not None else None,
                        "Day 2 p-value": round(d2_p, 4) if d2_p is not None else None,
                        "Day 7 log2FC": round(d7_fc, 2) if d7_fc is not None else None,
                        "Day 7 p-value": round(d7_p, 4) if d7_p is not None else None,
                        "Day 14 log2FC": round(d14_fc, 2) if d14_fc is not None else None,
                        "Day 14 p-value": round(d14_p, 4) if d14_p is not None else None,
                    })
    
    if results:
        display_df = pd.DataFrame(results)
        st.success(f"✅ Found {len(display_df)} matches")
        
        # Column order
        cols = ["Protein Accession", "Gene", "Protein Name"]
        if show_model: cols.append("Model")
        if show_species: cols.append("Species")
        if show_strain: cols.append("Strain")
        if show_gender: cols.append("Gender")
        if show_tissue: cols.append("Tissue")
        
        cols.extend(["Day 2 log2FC", "Day 2 p-value", "Day 7 log2FC", "Day 7 p-value", 
                    "Day 14 log2FC", "Day 14 p-value"])
        
        # Apply color styling
        def color_fc(val):
            if pd.isna(val):
                return ''
            color = 'background-color: #90EE90' if val > 0 else 'background-color: #FFB3B3'
            return color
        
        styled = display_df[cols].style.map(color_fc, subset=[c for c in cols if "log2FC" in c])
        
        st.dataframe(styled, use_container_width=True, hide_index=True)
        st.download_button("📥 Download Results", display_df.to_csv(index=False), "idea_results.csv")
    else:
        st.warning("No matching targets found.")

# Footer
st.markdown("---")
st.markdown("**© Ichor Life Sciences, Inc.** All rights reserved. This tool and its contents are proprietary to Ichor Life Sciences, Inc. Unauthorized use or distribution is prohibited. [www.ichorlifesciences.com](https://www.ichorlifesciences.com)")
