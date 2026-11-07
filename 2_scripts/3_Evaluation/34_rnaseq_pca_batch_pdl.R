
suppressPackageStartupMessages({
  library(DESeq2)
})

#-----------------------------------------------------------#
# Project paths
#-----------------------------------------------------------#

# Resolve project root relative to this R script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(
    sub("^--file=", "", file_arg[1])
  )
  script_dir <- dirname(script_path)
} else {
  stop("Could not determine script location. Run this script using Rscript.")
}

project_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  mustWork = TRUE
)

# Allow users to override input data location
base <- Sys.getenv(
  "TELOS_RNASEQ_BASE",
  unset = "/tmp/hayflick_review"
)

count_file <- file.path(
  base, "data", "RNAseq", "count_tables", "RS_raw_counts.txt"
)

meta_file <- file.path(
  base, "data_processing", "RNAseq", "DEseq", "RS_column_data.txt"
)

# Output directories
results_dir <- file.path(project_root, "3_results")

outdir <- file.path(results_dir, "figures")
tabdir <- file.path(results_dir, "tables", "rnaseq")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

# Check input files
if (!file.exists(count_file)) {
  stop("Count file not found: ", count_file)
}

if (!file.exists(meta_file)) {
  stop("Metadata file not found: ", meta_file)
}

#-----------------------------------------------------------#
# Read counts and metadata
#-----------------------------------------------------------#

counts <- read.delim(
  count_file,
  row.names = 1,
  check.names = FALSE
)

meta <- read.delim(
  meta_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(counts) <- gsub('"', '', colnames(counts))
rownames(counts) <- gsub('"', '', rownames(counts))
colnames(meta) <- gsub('"', '', colnames(meta))

meta$sample <- gsub('"', '', meta$sample)
meta$batch <- gsub('"', '', meta$batch)
meta$TP <- gsub('"', '', meta$TP)

meta$PDL <- as.numeric(
  sub(".*PDL([0-9]+).*", "\\1", meta$sample)
)

stopifnot(!anyDuplicated(meta$sample))
stopifnot(all(meta$sample %in% colnames(counts)))

# Align metadata to count columns
meta <- meta[match(colnames(counts), meta$sample), ]

stopifnot(all(meta$sample == colnames(counts)))

#-----------------------------------------------------------#
# Prepare count matrix
#-----------------------------------------------------------#

# Raw counts must be non-negative integers for DESeq2
count_matrix <- as.matrix(round(counts))
storage.mode(count_matrix) <- "integer"

if (any(count_matrix < 0)) {
  stop("Negative counts detected.")
}

rownames(meta) <- meta$sample
meta$batch <- factor(meta$batch)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = meta,
  design = ~ batch
)

# Filter very low-count genes before transformation
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

#-----------------------------------------------------------#
# Variance stabilizing transformation and PCA
#-----------------------------------------------------------#

# VST for visualisation only; this is not differential expression testing
vsd <- vst(dds, blind = TRUE)

pca <- prcomp(
  t(assay(vsd)),
  center = TRUE,
  scale. = FALSE
)

percent <- 100 * (
  pca$sdev^2 / sum(pca$sdev^2)
)

#-----------------------------------------------------------#
# Export gene loadings
#-----------------------------------------------------------#

loadings <- data.frame(
  gene = rownames(pca$rotation),
  PC1_loading = pca$rotation[, 1],
  PC2_loading = pca$rotation[, 2]
)

write.csv(
  loadings,
  file.path(tabdir, "RNAseq_PCA_gene_loadings.csv"),
  row.names = FALSE
)

#-----------------------------------------------------------#
# Export PCA coordinates
#-----------------------------------------------------------#

pca_df <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  batch = meta[rownames(pca$x), "batch"],
  PDL = meta[rownames(pca$x), "PDL"],
  TP = meta[rownames(pca$x), "TP"]
)

write.csv(
  pca_df,
  file.path(tabdir, "RNAseq_PCA_coordinates_batch_PDL.csv"),
  row.names = FALSE
)

#-----------------------------------------------------------#
# PCA plotting
#-----------------------------------------------------------#

# Base R plotting to avoid dependency on ggplot2
png_file <- file.path(outdir, "RNAseq_PCA_batch_PDL.png")
pdf_file <- file.path(outdir, "RNAseq_PCA_batch_PDL.pdf")

draw_pca <- function() {

  batch_colors <- c(
    "A" = "steelblue3",
    "B" = "darkorange2"
  )

  pch_values <- c(
    "A" = 16,
    "B" = 17
  )

  plot(
    pca_df$PC1,
    pca_df$PC2,
    col = batch_colors[as.character(pca_df$batch)],
    pch = pch_values[as.character(pca_df$batch)],
    cex = 1.3,
    xlab = sprintf("PC1 (%.1f%% variance)", percent[1]),
    ylab = sprintf("PC2 (%.1f%% variance)", percent[2]),
    main = "WI-38 RNA-seq PCA: batch and population doubling",
    bty = "l"
  )

  text(
    pca_df$PC1,
    pca_df$PC2,
    labels = paste0("PDL", pca_df$PDL),
    pos = 3,
    cex = 0.65,
    offset = 0.4
  )

  legend(
    "topright",
    legend = c("Batch A", "Batch B"),
    col = batch_colors[c("A", "B")],
    pch = pch_values[c("A", "B")],
    bty = "n"
  )

  grid(col = "grey85", lty = 3)
}

# Save PNG
png(
  png_file,
  width = 2400,
  height = 1800,
  res = 200
)

draw_pca()
dev.off()

# Save PDF
pdf(
  pdf_file,
  width = 10,
  height = 7.5
)

draw_pca()
dev.off()

#-----------------------------------------------------------#
# Summary
#-----------------------------------------------------------#

cat("PCA completed successfully.\n")
cat("Samples:", nrow(pca_df), "\n")
cat("Genes retained:", nrow(dds), "\n")

cat(sprintf(
  "PC1 variance: %.2f%%\n",
  percent[1]
))

cat(sprintf(
  "PC2 variance: %.2f%%\n",
  percent[2]
))

cat("PNG:", png_file, "\n")
cat("PDF:", pdf_file, "\n")

cat(
  "Coordinates:",
  file.path(tabdir, "RNAseq_PCA_coordinates_batch_PDL.csv"),
  "\n"
)

cat("\nSamples per batch:\n")
print(table(meta$batch))

cat("\nPDL by batch:\n")
print(table(meta$PDL, meta$batch))
