
    
import scipy.stats as st
import pandas as pd
import numpy as np
import warnings
warnings.simplefilter('ignore')

# --- Load Data ---
df_ogt = pd.read_csv("../data/csv/OGT.csv", index_col=0)
df_ogt = df_ogt[df_ogt["domain"] == "Bacteria"]

df_gene = pd.read_csv("../data/csv/gene_count.csv", index_col=0)


df_merged = df_gene.merge(
    df_ogt[["species", "phylum", "accession", "genus", "OGT"]],
    on="species",
    how="inner"
)




def search_gene(df_merged: pd.DataFrame, phylum: str):
   
   
    df = df_merged[df_merged["phylum"] != phylum].copy()
    ko_columns = [c for c in df_merged.columns if c.startswith("K")]
    

    gene_hot = df[df["OGT"] >= 20][ko_columns].copy()
    gene_cold = df[df["OGT"] < 20][ko_columns].copy()

    gene_hot = (gene_hot >= 1).astype(int)
    #print(gene_hot)
    gene_cold = (gene_cold >= 1).astype(int)

    

    results = []
    hot_sum = gene_hot.sum()
    #print(hot_sum)
    cold_sum = gene_cold.sum()

    # --- Fisher's exact test ---
    for ko in ko_columns:
        a = hot_sum[ko]
        b = gene_hot.shape[0] - a
        c = cold_sum[ko]
        d = gene_cold.shape[0] - c

        table = np.array([[a, b], [c, d]])
        pvalue = st.fisher_exact(table)[1]

        results.append({
            "ko": ko,
            "pvalue": pvalue,
            "data": f"{a},{b},{c},{d}"
        })

    df_pvalue = pd.DataFrame(results).sort_values("pvalue")
    df_pvalue.iloc[:50, 0:1].to_csv(f"../data/gene_list_without_phylum/method2/gene_{phylum}_cold.csv", index=False)


    return df_pvalue

for phylum in df_merged["phylum"].value_counts().keys()[:]:
    search_gene(df_merged,phylum)