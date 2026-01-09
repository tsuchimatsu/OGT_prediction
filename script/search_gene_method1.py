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

    hot_list = []
    cold_list = []


    for genus in df["genus"].value_counts().keys():
        sub = df[df["genus"] == genus]
        sub_sorted = sub.sort_values("OGT", ascending=False)

        hi, lo = 0, len(sub_sorted) - 1
        cnt = 0

        while lo > hi:
            t_hot = sub_sorted.iloc[hi]["OGT"]
            t_cold = sub_sorted.iloc[lo]["OGT"]

            if (t_hot - t_cold < 10) or (cnt >= 3):
                break

            hot_list.append(sub_sorted.iloc[hi]["accession"])
            cold_list.append(sub_sorted.iloc[lo]["accession"])

            hi += 1
            lo -= 1
            cnt += 1

  
    ko_columns = [c for c in df_merged.columns if c.startswith("K")]
    

    gene_hot = df[df["accession"].isin(hot_list)][ko_columns].copy()
    
    gene_cold = df[df["accession"].isin(cold_list)][ko_columns].copy()

   
    gene_hot = (gene_hot >= 1).astype(int)
    #print(gene_hot)
    gene_cold = (gene_cold >= 1).astype(int)

    num = gene_hot.shape[0]

    results = []
    hot_sum = gene_hot.sum()
    #print(hot_sum)
    cold_sum = gene_cold.sum()

    # --- Fisher's exact test ---
    for ko in ko_columns:
        a = hot_sum[ko]
        b = num - a
        c = cold_sum[ko]
        d = num - c

        table = np.array([[a, b], [c, d]])
        pvalue = st.fisher_exact(table)[1]

        results.append({
            "ko": ko,
            "pvalue": pvalue,
            "data": f"{a},{b},{c},{d}"
        })

    df_pvalue = pd.DataFrame(results).sort_values("pvalue")
    df_pvalue.iloc[:50, 0:1].to_csv(f"../data/gene_list_without_phylum/method1/gene_{phylum}.csv", index=False)


    return df_pvalue

for phylum in df_merged["phylum"].value_counts().keys()[:]:
    search_gene(df_merged,phylum)