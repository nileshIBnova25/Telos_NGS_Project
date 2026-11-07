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
if (!requireNamespace("limma", quietly = TRUE)) stop("Install limma: BiocManager::install('limma')")
library(ggplot2)
library(ggrepel)
library(limma)

# Prepare unique gene rows and numeric expression matrix.
expr_df <- tpm |>
  dplyr::select(gene, dplyr::all_of(sample_cols)) |>
  dplyr::filter(!is.na(gene), as.character(gene) != "") |>
  dplyr::distinct(gene, .keep_all = TRUE)

gene_names <- make.unique(as.character(expr_df$gene))
mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
storage.mode(mat) <- "numeric"
rownames(mat) <- gene_names

if (anyNA(mat)) warning("Expression matrix contains missing values; rows with missing values may be excluded.")
keep <- rowMeans(mat, na.rm = TRUE) > 1 & rowSums(is.finite(mat)) == ncol(mat)
mat <- mat[keep, , drop = FALSE]
if (nrow(mat) < 2) stop("Too few genes remain after expression filtering.")

logmat <- log2(mat + 0.001)
group <- factor(c(rep("Young", length(young)), rep("RS", length(rs))),
                levels = c("Young", "RS"))
design <- model.matrix(~ group)
fit <- limma::lmFit(logmat, design)
fit <- limma::eBayes(fit)

res <- limma::topTable(fit, coef = "groupRS", number = Inf, sort.by = "none") |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    direction = dplyr::case_when(
      adj.P.Val < 0.05 & logFC >= 1 ~ "Up in RS",
      adj.P.Val < 0.05 & logFC <= -1 ~ "Down in RS",
      TRUE ~ "Not significant"
    ),
    neglog10_p = -log10(pmax(P.Value, 1e-300))
  )

write.csv(res, file.path(out_dir, "exploratory_limma_results.csv"), row.names = FALSE)

label_genes <- res |>
  dplyr::filter(adj.P.Val < 0.05) |>
  dplyr::arrange(adj.P.Val) |>
  dplyr::slice_head(n = 12)

p <- ggplot(res, aes(x = logFC, y = neglog10_p, color = direction)) +
  geom_point(alpha = 0.7, size = 1.25) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  ggrepel::geom_text_repel(data = label_genes, aes(label = gene),
                           size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c("Up in RS" = "#D73027",
                                "Down in RS" = "#4575B4",
                                "Not significant" = "grey75")) +
  labs(
    title = "Young vs late/RS WI-38 cells",
    subtitle = "Exploratory limma analysis of log2(TPM + 0.001); not a count-based DE analysis",
    x = "Log2 fold change (RS vs Young)",
    y = "-Log10 nominal p-value",
    color = "Classification"
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "volcano_plot.png"), p, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "volcano_plot.pdf"), p, width = 8, height = 6)
print(p)
cat("Saved volcano plot and exploratory results to:", out_dir, "\n")
cat("NOTE: TPM-based p-values are exploratory, not definitive differential-expression evidence.\n")
