import streamlit as st
import pandas as pd
from pathlib import Path

st.set_page_config(page_title="Ichor Life Sciences • IDEA", layout="wide")

st.markdown("<h1 style='color:#1a3c6e; text-align:center;'>Ichor Life Sciences</h1>", unsafe_allow_html=True)
st.subheader("Differential Expression Atlas (IDEA)")

# Load all data-*.xlsx files
data_files = list(Path(".").glob("data-*.xlsx"))
st.success(f"✅ Found {len(data_files)} model dataset(s)")

models = {}
for file in data_files:
    model_name = file.stem.replace("data-", "")
    xls = pd.ExcelFile(file)
    
    # Cover sheet (key-value format)
    cover = pd.read_excel(xls, "cover", header=None)
    meta = dict(zip(cover.iloc[:, 0].astype(str).str.strip(), cover.iloc[:, 1].astype(str).str.strip()))
    
    log2 = pd.read_excel(xls, "log2")
    
    models[model_name] = {
        "meta": meta,
        "log2": log2
    }

# Search
query = st.text_input("🔍 Search by Protein Accession, Gene, or Protein Name", 
                     placeholder="Ca1, Alb, Gapdh, Col1a1")

if query and models:
    terms = [t.strip().upper() for t in query.split(",") if t.strip()]
    results = []
    
    for model_name, data in models.items():
        df = data["log2"]
        meta = data["meta"]
        
        # Search across Accession (col 0), Gene (col 1), Protein Name (col 2)
        mask = (df.iloc[:, 0].astype(str).str.upper().str.contains('|'.join(terms), na=False)) | \
               (df.iloc[:, 1].astype(str).str.upper().str.contains('|'.join(terms), na=False)) | \
               (df.iloc[:, 2].astype(str).str.upper().str.contains('|'.join(terms), na=False))
        
        hits = df[mask].copy()
        if len(hits) > 0:
            for _, row in hits.iterrows():
                results.append({
                    "Protein Accession": row.iloc[0],
                    "Gene": row.iloc[1],
                    "Protein Name": row.iloc[2],
                    "Model": model_name,
                    "Species": meta.get("Species", ""),
                    "Strain": meta.get("Strain", ""),
                    "Gender": meta.get("Gender", ""),
                    "Tissue": meta.get("Tissue", ""),
                    "Day 2 vs NAIVE (log2FC)": round(row.iloc[3], 2) if len(row) > 3 else None,
                    "Day 7 vs NAIVE (log2FC)": round(row.iloc[4], 2) if len(row) > 4 else None,
                    "Day 14 vs NAIVE (log2FC)": round(row.iloc[5], 2) if len(row) > 5 else None,
                })
    
    if results:
        display_df = pd.DataFrame(results)
        st.success(f"✅ Found {len(display_df)} matches")
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        st.download_button("📥 Download Results", display_df.to_csv(index=False), "idea_results.csv")
    else:
        st.warning("No matching targets found.")

st.caption("Multi-model support active. Add more `data-*.xlsx` files as needed.")
