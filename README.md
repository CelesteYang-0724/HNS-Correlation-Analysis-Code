# Figure 6A–B: Bacteria–Metabolite Correlation Analysis

This repository provides the R script used for bacteria–metabolite correlation analysis and visualization in Figures 6A and 6B.

## Analysis workflow

- Spearman correlation analysis between bacterial genera and metabolites.
- Hierarchical clustering heatmap of bacteria–metabolite correlations (Figure 6A).
- Sankey visualization of significant bacteria–metabolite associations (Figure 6B).

## Software

- R version 4.5.2
- Required R packages are specified in the analysis script.

## Repository structure

```text
HNS-Correlation-Analysis-Code/
├── Figure6_Correlation_Analysis.R
└── README.md
```

## Data availability

The analysis requires three processed input datasets:

- `属水平肠道菌群.csv`: Bacterial genus-level abundance matrix.
- `HNS代谢物-neg.csv`: Negative-ion metabolomics intensity matrix.
- `HNS代谢物-pos.csv`: Positive-ion metabolomics intensity matrix.

These datasets are not included in this public code repository. The corresponding archived release (v1.0.0), including code and processed datasets, is registered in Zenodo (https://doi.org/10.5281/zenodo.24475252). Access to the archived files is currently restricted pending publication.

## Usage

Place the required CSV files in the `data/` directory, verify the input file paths in the R script, and run:

```r
source("Figure6_Correlation_Analysis.R")
```

The code is publicly available; reproduction of the analyses requires access to the corresponding input datasets.
