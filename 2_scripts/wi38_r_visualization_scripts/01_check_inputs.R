#!/usr/bin/env Rscript
# Shared setup is repeated in each standalone script so each can run independently.
options(stringsAsFactors = FALSE)

project_root <- getwd()
gt_path <- file.path(project_root, "1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv")
tpm_path <- file.path(project_root, "1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx")
out_dir <- file.path(project_root, "1_data_raw/4_rnaseq/outputs/visualizations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

young <- c("PDL20_TP1_A", "PDL20_TP1_B", "PDL20_TP1_C")
rs <- c("PDL46_TP6_A", "PDL46_TP6_B", "PDL50_TP7_A", "PDL50_TP7_B")
sample_cols <- c(young, rs)

if (!file.exists(gt_path)) stop("Ground-truth CSV not found: ", gt_path, "\nRun Rscript from the project root or edit gt_path.")
if (!file.exists(tpm_path)) stop("TPM workbook not found: ", tpm_path, "\nEdit tpm_path in the script if your workbook is elsewhere.")

if (!requireNamespace("readxl", quietly = TRUE)) stop("Missing package readxl. Install with install.packages('readxl').")
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Missing package dplyr. Install with install.packages('dplyr').")
if (!requireNamespace("tibble", quietly = TRUE)) stop("Missing package tibble. Install with install.packages('tibble').")

library(readxl)
library(dplyr)
library(tibble)

gt <- read.csv(gt_path, check.names = FALSE)
tpm <- readxl::read_excel(tpm_path, sheet = "RS_TC_raw_TPM")

# Identify likely gene identifier and log2FC columns.
gene_candidates <- c("gene_id_or_symbol", "gene", "gene_symbol", "gene_name", "Gene", "GeneSymbol", "gene_id", "Gene ID")
fc_candidates <- c("log2FC_rs_vs_young", "log2FC", "log2FoldChange", "logFC", "log2_fc_rs_vs_young")
gene_col <- gene_candidates[gene_candidates %in% names(gt)][1]
fc_col <- fc_candidates[fc_candidates %in% names(gt)][1]

if (is.na(gene_col) || is.na(fc_col)) {
  message("Ground-truth columns detected: ", paste(names(gt), collapse = ", "))
  stop("Could not automatically identify gene and log2FC columns. Edit gene_candidates/fc_candidates in this script to match names(gt).")
}
if (!"gene" %in% names(tpm)) {
  message("TPM workbook columns detected: ", paste(names(tpm), collapse = ", "))
  stop("Expected a gene column named 'gene' in the TPM sheet. Edit the scripts to match your workbook.")
}
missing_samples <- setdiff(sample_cols, names(tpm))
if (length(missing_samples)) {
  message("TPM workbook columns detected: ", paste(names(tpm), collapse = ", "))
  stop("Missing expected sample columns: ", paste(missing_samples, collapse = ", "),
       "\nEdit young/rs sample names in the script to match your workbook.")
}

cat("\n--- Ground-truth table ---\n")
cat("Path:", gt_path, "\n")
cat("Dimensions:", nrow(gt), "rows x", ncol(gt), "columns\n")
cat("Columns:\n"); print(names(gt))
cat("\nFirst rows:\n"); print(utils::head(gt, 6))
cat("\nDetected gene column:", gene_col, "\nDetected log2FC column:", fc_col, "\n")

cat("\n--- TPM workbook ---\n")
cat("Path:", tpm_path, "\n")
cat("Sheet: RS_TC_raw_TPM\n")
cat("Dimensions:", nrow(tpm), "rows x", ncol(tpm), "columns\n")
cat("Columns:\n"); print(names(tpm))
cat("\nSample columns found:\n"); print(sample_cols)

cat("\n--- Output directory ---\n", out_dir, "\n", sep="")
cat("\nInput check complete. If the expected sample names or column names differ, edit the scripts before continuing.\n")
