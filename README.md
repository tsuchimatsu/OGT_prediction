# OGT Prediction from Genomic Features

This repository contains scripts and data for calculating genomic features and predicting the optimum growth temperature (OGT) of bacteria.

## Environment Setup

To create the conda environment, run:

```bash
conda env create -f env/environment.yml
```

## Data Description

CSV files calculated in this study are located in the `/data/csv/` directory:

- `calculated_features_archaea.csv.zip`: Genomic features of archaea.
- `calculated_features_bacteria.csv.zip`: Genomic features of bacteria.
- `gene_count.csv.zip`: The Number of each gene for each bacterial species.
- `OGT.csv`: Optimum growth temperature data for both bacteria and archaea.
- `bac_asr.csv`: Ancestral OGT of most recent common ancestor of the genus (can be created by calc_ancestral_OGT.ipynb)

Gene lists used for regression are located in the `/data/gene_list_without_phylum` directory:

Tree files used in this study are located in the `/data/tree/` directory:

- `ar53_r207_selected.tree`: tree of archaea.
- `bac120_r207_selected.tree`: tree of bacteria.

## Scripts

- `/script/calc_features.py`: Calculates genomic features from input genomes (example genome data available in `/data/genome_example`).
- `/script/predict_OGT.ipynb`: Jupyter notebook for predicting OGT using genomic features, with or without gene presence/absence information.
- `/script/search_genes_method1.py`: Searches genes using method1
- `/script/search_genes_method2.py`: Searches genes using method2
- `/script/calc_ancestral_OGT.ipynb`: Calculates ancestral OGT using castor (please use R 4.2 and install required packages)

