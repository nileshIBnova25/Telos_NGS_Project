# RNA-seq workflow: GSE175533

This folder contains the exploratory RNA-seq workflow for comparing young/proliferating WI-38 fibroblasts with late-passage replicative-senescence (RS) cells.

## Scripts

1. `1_rna_inspect_tpm_table.py` — inspects Excel sheet names and exact column headers, writes local inspection reports, and creates `1_data_raw/4_rnaseq/rna_config.json` if it does not already exist.
2. `2_rna_make_groundtruth.py` — after you manually confirm the exact gene column and replicate columns in that JSON, computes group mean TPM, `log2(TPM + 1)`, and `log2FC_rs_vs_young`.
3. `3_rna_qc_summary.py` — writes basic QC statistics plus top positive/negative TPM log2-fold-change lists.

`rna_utils.py` contains shared path, logging, and completion-marker helpers.

## Run order

Run commands from the repository root (`Telos_NGS_Project`). Use the Python environment in which `pandas`, `numpy`, and `openpyxl` are installed.

```bash
python 2_scripts/2_RNA/1_rna_inspect_tpm_table.py
```

Inspect the log and these local reports:

- `1_data_raw/4_rnaseq/inspection/workbook_sheets.csv`
- `1_data_raw/4_rnaseq/inspection/columns_by_sheet.csv`
- `1_data_raw/4_rnaseq/inspection/sample_column_matches.csv`

Then edit `1_data_raw/4_rnaseq/rna_config.json` and set:

- `gene_column`: exact gene identifier column in `RS_TC_raw_TPM`
- `young_columns`: exact WI-38 PDL20 replicate columns (expected A/B/C if those are present)
- `rs_columns`: exact WI-38 PDL46/PDL50 replicate columns (expected A/B for each if those are present)

Do not include `hTERT`, `RIS`, or `CD` columns in this primary RS comparison. Do not use `RS_TC_batch_corrected_TPM` for this calculation; the project notes reserve that sheet for visualization.

```bash
python 2_scripts/2_RNA/2_rna_make_groundtruth.py
python 2_scripts/2_RNA/3_rna_qc_summary.py
```

## Outputs (local only)

- `1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv`
- `1_data_raw/4_rnaseq/qc/rnaseq_qc_summary.csv`
- `1_data_raw/4_rnaseq/qc/top_up_genes.csv`
- `1_data_raw/4_rnaseq/qc/top_down_genes.csv`
- Logs and `R1`/`R2`/`R3` `.done.txt` markers under `1_data_raw/4_rnaseq/`

## Interpretation caveats

- TPM mean and log2 fold change are descriptive summaries. They do not replace an appropriately designed statistical differential-expression analysis on raw counts.
- The top-up/top-down files are ranked by TPM log2 fold change and are not significance-filtered.
- Check sample metadata before running; these scripts deliberately do not guess sample columns from keyword matches.
- If the gene identifier column contains Ensembl IDs rather than gene symbols, map IDs to symbols before interpreting named genes such as `YAP1`, `TEAD1`, or `TGFB2`.
