import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path

VERSION = "1.2"

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

# Logo
logo_path = "logo.png"
if Path(logo_path).exists():
    st.image(logo_path, width=320)
else:
    st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)

st.subheader("Ichor Differential Expression Atlas (IDEA)")
st.caption(f"**Pre-Clinical Model Explorer** — Version {VERSION}")

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
fc_filter = st.sidebar.slider("|log2FC| ≥ (0 = show all)", 0.0, 5.0, 0.0, 0.1)
p_filter = st.sidebar.slider("p-value < (0.1 = show all)", 0.0001, 0.1, 0.1, 0.001)

st.sidebar.header("Display Columns")
show_model = st.sidebar.checkbox("Model", value=True)
show_gene = st.sidebar.checkbox("Gene", value=True)
show_protein = st.sidebar.checkbox("Protein Name", value=True)
show_species = st.sidebar.checkbox("Species", value=True)
show_strain = st.sidebar.checkbox("Strain", value=True)
show_gender = st.sidebar.checkbox("Gender", value=True)
show_tissue = st.sidebar.checkbox("Tissue", value=True)
show_pvalues = st.sidebar.checkbox("P-values", value=False)

st.sidebar.header("Color Coding")
fc_color_thresh = st.sidebar.slider("Color |log2FC| threshold", 0.0, 5.0, 1.0, 0.1)

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
                d2_fc = round(float(row.iloc[3]), 2) if len(row) > 3 and pd.notna(row.iloc[3]) else None
                d7_fc = round(float(row.iloc[4]), 2) if len(row) > 4 and pd.notna(row.iloc[4]) else None
                d14_fc = round(float(row.iloc[5]), 2) if len(row) > 5 and pd.notna(row.iloc[5]) else None
                
                d2_p = round(float(p_df.iloc[idx, 3]), 4) if len(p_df.columns) > 3 and pd.notna(p_df.iloc[idx, 3]) else None
                d7_p = round(float(p_df.iloc[idx, 4]), 4) if len(p_df.columns) > 4 and pd.notna(p_df.iloc[idx, 4]) else None
                d14_p = round(float(p_df.iloc[idx, 5]), 4) if len(p_df.columns) > 5 and pd.notna(p_df.iloc[idx, 5]) else None
                
                include = True
                if fc_filter > 0:
                    include &= any(abs(x or 0) >= fc_filter for x in [d2_fc, d7_fc, d14_fc])
                if p_filter < 0.1:
                    include &= any((x or 1) < p_filter for x in [d2_p, d7_p, d14_p])
                
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
                        "Day 2 log2FC": d2_fc,
                        "Day 2 p-value": d2_p,
                        "Day 7 log2FC": d7_fc,
                        "Day 7 p-value": d7_p,
                        "Day 14 log2FC": d14_fc,
                        "Day 14 p-value": d14_p,
                    })
    
    if results:
        display_df = pd.DataFrame(results)
        st.success(f"✅ Found {len(display_df)} matches")
        
        cols = ["Protein Accession", "Gene", "Protein Name"]
        if show_model: cols.append("Model")
        if show_species: cols.append("Species")
        if show_strain: cols.append("Strain")
        if show_gender: cols.append("Gender")
        if show_tissue: cols.append("Tissue")
        
        fc_p_cols = ["Day 2 log2FC", "Day 2 p-value", "Day 7 log2FC", "Day 7 p-value", 
                    "Day 14 log2FC", "Day 14 p-value"]
        
        if show_pvalues:
            cols.extend(fc_p_cols)
        else:
            cols.extend([c for c in fc_p_cols if "log2FC" in c])
        
        # Color coding
        def color_fc(val):
            if pd.isna(val):
                return ''
            if abs(val) >= fc_color_thresh:
                return 'background-color: #90EE90' if val > 0 else 'background-color: #FFB3B3'
            return ''
        
        styled = display_df[cols].style.map(color_fc, subset=[c for c in cols if "log2FC" in c])
        
        st.dataframe(styled, use_container_width=True, hide_index=True)
        
        # Volcano Plot
        if st.button("Show Volcano Plot (Day 7 vs NAIVE)"):
            plot_df = display_df.copy()
            plot_df['-log10(p-value)'] = -pd.np.log10(plot_df['Day 7 p-value'].replace(0, 1e-10))
            fig = px.scatter(plot_df, x="Day 7 log2FC", y="-log10(p-value)", 
                            hover_data=["Gene", "Protein Name", "Model"],
                            title="Volcano Plot - Day 7 vs NAIVE",
                            labels={"Day 7 log2FC": "log2 Fold Change"})
            fig.add_hline(y=-pd.np.log10(0.05), line_dash="dash", line_color="gray")
            fig.add_vline(x=1, line_dash="dash", line_color="gray")
            fig.add_vline(x=-1, line_dash="dash", line_color="gray")
            st.plotly_chart(fig, use_container_width=True)
        
        st.download_button("📥 Download Results", display_df.to_csv(index=False), "idea_results.csv")
    else:
        st.warning("No matching targets found.")

# Footer
st.markdown("---")
st.markdown("**© Ichor Life Sciences, Inc.** All rights reserved. This tool and its contents are proprietary to Ichor Life Sciences, Inc. Unauthorized use or distribution is prohibited. [www.ichorlifesciences.com](https://www.ichorlifesciences.com)")
