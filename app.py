import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)
st.subheader("Differential Expression Atlas (IDEA)")

# Load all data-* files
data_files = list(Path(".").glob("data-*.xlsx"))
st.success(f"✅ Found {len(data_files)} model dataset(s)")

models = {}
for file in data_files:
    model_name = file.stem.replace("data-", "")
    xls = pd.ExcelFile(file)
    
    cover = pd.read_excel(xls, "cover")
    log2 = pd.read_excel(xls, "log2")
    pval = pd.read_excel(xls, "p")
    
    models[model_name] = {
        "cover": cover,
        "log2": log2,
        "p": pval
    }

# Search
query = st.text_input("🔍 Search by Protein Accession, Gene, or Protein Name", 
                     placeholder="Ca1, Alb, Gapdh, P47739")

if query and models:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    results = []
    
    for model_name, data in models.items():
        df = data["log2"]
        cover = data["cover"]
        
        mask = (df.iloc[:, 0].astype(str).str.upper().str.contains('|'.join(terms), na=False)) | \
               (df.iloc[:, 1].astype(str).str.upper().str.contains('|'.join(terms), na=False)) | \
               (df.iloc[:, 2].astype(str).str.upper().str.contains('|'.join(terms), na=False))
        
        hits = df[mask].copy()
        if len(hits) > 0:
            meta = cover.iloc[0].to_dict()
            for _, row in hits.iterrows():
                results.append({
                    "Model": model_name,
                    "Protein Accession": row.iloc[0],
                    "Gene": row.iloc[1],
                    "Protein Name": row.iloc[2],
                    "Species": meta.get("Species", ""),
                    "Strain": meta.get("Strain", ""),
                    "Gender": meta.get("Gender", ""),
                    "Tissue": meta.get("Tissue", ""),
                    "Day 2 log2FC": row.iloc[3],
                    "Day 7 log2FC": row.iloc[4],
                    "Day 14 log2FC": row.iloc[5],
                })
    
    if results:
        display_df = pd.DataFrame(results)
        st.success(f"✅ Found {len(display_df)} matches across {len(models)} model(s)")
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.download_button("📥 Download Results", display_df.to_csv(index=False), "idea_results.csv")
    else:
        st.warning("No matching targets found.")

st.caption("Scalable multi-model support active. Add more `data-*.xlsx` files to expand.")
