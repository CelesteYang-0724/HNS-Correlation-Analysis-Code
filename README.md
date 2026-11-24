# Gut microbiota correlation analysis — R code

This public repository contains reproducible **analysis code only**, for two separate workflows in the manuscript. The input data are **not included** in this public code repository.

## 1. Bacteria–metabolite correlation (Figures 6A–B)

**Script:** `Figure6_Correlation_Analysis.R`

- Spearman correlation analysis of genus-level abundances and positive-/negative-ion metabolomics features.
- Hierarchical clustering heatmap (Figure 6A) and Sankey/alluvial visualization (Figure 6B).
- Developed with R 4.5.2. See the script for required R packages and thresholds.
- Input files required locally under `data/`: `属水平肠道菌群.csv`, `HNS代谢物-neg.csv`, `HNS代谢物-pos.csv`.

```r
source("Figure6_Correlation_Analysis.R", encoding = "UTF-8")
```

The corresponding **earlier v1.0.0 release** is permanently recorded at [Zenodo DOI: 10.5281/zenodo.24475252](https://doi.org/10.5281/zenodo.24475252). The archived files are **currently restricted**. This DOI is for the original Figure 6A–B release, not for the new phenotype correlation script.

## 2. Bacteria–biochemical/behavioral phenotype correlation

**Script:** `Microbiome_Phenotype_Correlation.R`

- Spearman correlation analysis of bacterial genera against inflammatory, biochemical and behavioral measurements.
- Saves PDF and 600-dpi PNG heatmaps and both complete and significant association tables.
- Uses base R, tested target environment R 4.5.2.
- Input files required locally under `data/`: `属水平肠道菌群n=6.csv` and `炎症因子及行为学指标.csv`.

```r
source("Microbiome_Phenotype_Correlation.R", encoding = "UTF-8")
```

The phenotype analysis **matches sample IDs between files**. It will not pair unmatched samples by row position. See [detailed instructions](README_Microbiome_Phenotype.md) before running or reusing results.

## Repository layout

```text
HNS-Correlation-Analysis-Code/
├── Figure6_Correlation_Analysis.R
├── Microbiome_Phenotype_Correlation.R
├── README.md
├── README_Microbiome_Phenotype.md
└── .gitignore
```

Users with authorized access to the inputs should create a **local** `data/` directory; generated files appear under `results/`. Both directories are excluded from Git tracking by `.gitignore`.

## Data accessibility and reproducibility

The **source/processed input datasets are not hosted in this public code repository**. Availability is governed by the study's data-sharing permissions and journal requirements. The archived Figure 6A–B dataset at the Zenodo DOI above currently has restricted file access; it does not contain the phenotype analysis datasets. Accordingly, code availability **does not imply that the analyses can be reproduced without obtaining the appropriate inputs**.
