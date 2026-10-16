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

if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Install pheatmap: install.packages('pheatmap')")

fc <- suppressWarnings(as.numeric(gt[[fc_col]]))
gene_ids <- as.character(gt[[gene_col]])
gt2 <- data.frame(gene_id = gene_ids, log2FC = fc) |>
  dplyr::filter(!is.na(gene_id), gene_id != "", is.finite(log2FC)) |>
  dplyr::arrange(dplyr::desc(abs(log2FC))) |>
  dplyr::distinct(gene_id, .keep_all = TRUE) |>
  dplyr::slice_head(n = 100)

# Match exact gene identifiers to TPM gene column.
tpm_gene <- as.character(tpm$gene)
idx <- match(gt2$gene_id, tpm_gene)
matched <- !is.na(idx)
if (!any(matched)) {
  stop("No ground-truth gene IDs matched tpm$gene. Check whether the two files use symbols vs Ensembl IDs. Ground-truth gene column: ",
       gene_col, "; TPM gene column: gene.")
}
if (sum(matched) < 10) warning("Only ", sum(matched), " genes matched. Verify identifier types.")

heat_df <- tpm[idx[matched], sample_cols, drop = FALSE]
heat_mat <- as.matrix(heat_df)
storage.mode(heat_mat) <- "numeric"
rownames(heat_mat) <- make.unique(gt2$gene_id[matched])

heat_log <- log2(heat_mat + 0.001)
heat_z <- t(scale(t(heat_log)))
heat_z <- heat_z[apply(heat_z, 1, function(x) all(is.finite(x))), , drop = FALSE]
if (nrow(heat_z) < 2) stop("Too few variable matched genes for heatmap.")

annotation <- data.frame(
  Condition = factor(c(rep("Young", length(young)), rep("RS", length(rs))),
                     levels = c("Young", "RS"))
)
rownames(annotation) <- sample_cols

png(file.path(out_dir, "heatmap_top100_genes.png"), width = 1800, height = 2200, res = 200)
pheatmap::pheatmap(
  heat_z,
  annotation_col = annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = nrow(heat_z) <= 60,
  fontsize = 8,
  main = "Top changing WI-38 genes\nRow Z-scores of log2(TPM + 0.001)",
  silent = TRUE
)
dev.off()

pdf(file.path(out_dir, "heatmap_top100_genes.pdf"), width = 9, height = 11)
pheatmap::pheatmap(
  heat_z,
  annotation_col = annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = nrow(heat_z) <= 60,
  fontsize = 8,
  main = "Top changing WI-38 genes\nRow Z-scores of log2(TPM + 0.001)"
)
dev.off()

cat("Matched", sum(matched), "of", nrow(gt2), "selected genes.\n")
cat("Saved heatmap to:", out_dir, "\n")
