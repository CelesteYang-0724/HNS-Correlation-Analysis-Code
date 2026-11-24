# Microbiome–phenotype correlation analysis

This directory contains an R script for **Spearman correlation analysis between bacterial genera and biochemical/behavioral phenotypes**. It is derived from the project's earlier script `菌群生化指标相关性分析_英文alpha_beta_600dpi(1).R` and has been adapted for public sharing. No source data are included.

## Required inputs (not publicly distributed)

Place these **two** CSV files in a local `data/` directory:

| File | Contents |
| --- | --- |
| `属水平肠道菌群n=6.csv` | Bacterial genus abundance matrix |
| `炎症因子及行为学指标.csv` | Biochemical, inflammatory and behavioral phenotype measurements |

Each CSV must contain **sample IDs in the first column**, with samples in rows and numeric variables in subsequent columns. Sample IDs must be unique, and the same sample must use the same ID in both files. The script pairs records **by sample ID**, not by their row positions. At least three matching samples are needed; check any warning about unmatched samples. If the datasets cannot be reliably paired, the correlation results are not valid.

## Installation and run

The script uses **base R** (developed for R 4.5.2), without additional CRAN packages.

From the repository root, create the `data/` directory locally and add the two CSV files. Then run:

```r
source("Microbiome_Phenotype_Correlation.R", encoding = "UTF-8")
```

If the files were saved in GBK encoding, change `CSV_ENCODING <- ""` near the top of the script to `CSV_ENCODING <- "GBK"`.

Outputs are written automatically to `results/`; the script **does not open figures in the RStudio plotting pane**.

## Outputs

| File | Description |
| --- | --- |
| `results/correlation_heatmap.pdf` | Vector-format heatmap |
| `results/correlation_heatmap_600dpi.png` | 600-dpi PNG heatmap |
| `results/all_correlations.csv` | All tested pairs: Spearman rho, P value, BH-adjusted P, sample count |
| `results/significant_correlations.csv` | Pairs with **uncorrected** P < 0.05, as in the original script |
| `results/correlation_matrix.csv` | Full correlation coefficient matrix |
| `results/p_value_matrix.csv` | Full unadjusted P-value matrix |
| `results/sessionInfo.txt` | R session and package information |

The **heatmap stars and significant-results filter use unadjusted P values** (`*` P < 0.05, `**` P < 0.01, `***` P < 0.001). The additional `P_adjust_BH` column is provided for multiple-testing transparency; it is **not** used to select the plotted stars. The heatmap color range saturates at Spearman rho ±0.6, matching the original script's display settings. A correlation does not establish causation.

## Important difference from the original local script

The original script assigned artificial IDs to duplicates and truncated tables to equal row counts before correlating. Those operations could pair different biological samples. The public version **rejects duplicate sample IDs** and matches records by their real identifiers; therefore, if the original CSVs were not correctly aligned, numerical results may differ. Verify the sample mapping and inspect the full outputs before using this version for manuscript figures or replacing earlier statistical results.

## Data availability

**Only code is published in the public GitHub repository.** The input CSVs are not distributed here. The separate Figure 6A–B bacteria–metabolite analysis archive, [Zenodo DOI: 10.5281/zenodo.24475252](https://doi.org/10.5281/zenodo.24475252), is currently restricted and **does not include this script's two phenotype input CSVs**. Do not treat that Zenodo record as an archive of this analysis.

Any future sharing of the inputs should comply with the study's authorization and the journal's data-availability requirements.
