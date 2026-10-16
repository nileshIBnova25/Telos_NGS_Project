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

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2: install.packages('ggplot2')")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install ggrepel: install.packages('ggrepel')")
library(ggplot2)
library(ggrepel)

expr_df <- tpm |>
  dplyr::select(gene, dplyr::all_of(sample_cols)) |>
  dplyr::filter(!is.na(gene), as.character(gene) != "") |>
  dplyr::distinct(gene, .keep_all = TRUE)

mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- make.unique(as.character(expr_df$gene))

keep <- rowMeans(mat, na.rm = TRUE) > 1 & rowSums(is.finite(mat)) == ncol(mat)
mat <- mat[keep, , drop = FALSE]
logmat <- log2(mat + 0.001)

gene_var <- apply(logmat, 1, var, na.rm = TRUE)
gene_var <- gene_var[is.finite(gene_var)]
if (length(gene_var) < 3) stop("Too few usable genes for PCA.")
top_n <- min(1000, length(gene_var))
top_var <- names(sort(gene_var, decreasing = TRUE))[seq_len(top_n)]

pca <- prcomp(t(logmat[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
  tibble::rownames_to_column("Sample") |>
  dplyr::mutate(
    Condition = factor(ifelse(Sample %in% young, "Young", "RS"),
                       levels = c("Young", "RS"))
  )

p <- ggplot(pca_df, aes(PC1, PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
  labs(
    title = "PCA of WI-38 RNA expression",
    subtitle = paste0("Top ", top_n, " most variable genes; log2(TPM + 0.001)"),
    x = sprintf("PC1 (%.1f%% variance)", variance[1]),
    y = sprintf("PC2 (%.1f%% variance)", variance[2])
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "pca_young_vs_rs.png"), p, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "pca_young_vs_rs.pdf"), p, width = 8, height = 6)
print(p)
cat("Saved PCA plots to:", out_dir, "\n")
