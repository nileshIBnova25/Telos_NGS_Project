# WI-38 RNA-seq visualization scripts

Run these scripts from your project root (the directory containing `1_data_raw/`).

## Files
- `01_check_inputs.R` — locates and inspects the input files/columns.
- `02_volcano_exploratory.R` — exploratory volcano plot from log2(TPM + pseudocount) using limma.
- `03_heatmap_top_genes.R` — heatmap of top genes by absolute ground-truth log2 fold change.
- `04_pca_samples.R` — PCA of the individual sample TPM profiles.
- `run_all.sh` — runs the scripts in order using `Rscript`.

## Before running
These scripts assume:
- Ground truth CSV: `1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv`
- TPM workbook: `1_data_raw/4_rnaseq/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx`
- TPM sheet: `RS_TC_raw_TPM`
- Sample columns: `PDL20_TP1_A`, `PDL20_TP1_B`, `PDL20_TP1_C`,
  `PDL46_TP6_A`, `PDL46_TP6_B`, `PDL50_TP7_A`, `PDL50_TP7_B`

The scripts try to detect the gene ID column and the log2FC column from common names. If they cannot, edit the variables noted in the script.

## Install packages
Run once:
```bash
Rscript -e 'install.packages(c("readxl","dplyr","ggplot2","pheatmap","ggrepel","tibble"))'
Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("limma")'
```

## Run
```bash
chmod +x run_all.sh
./run_all.sh
```

Or individually:
```bash
Rscript 01_check_inputs.R
Rscript 02_volcano_exploratory.R
Rscript 03_heatmap_top_genes.R
Rscript 04_pca_samples.R
```

Outputs are written to `1_data_raw/4_rnaseq/outputs/visualizations/`.

## Important scientific note
The volcano script is explicitly exploratory because it uses TPM values transformed to log2(TPM + 0.001), not raw counts or transcript-level quantifications. For inferential RNA-seq differential expression, use the appropriate raw counts/quantifications and a method such as DESeq2, edgeR, or limma-voom, with a design that reflects the actual experimental structure. Do not treat the exploratory p-values as definitive.
