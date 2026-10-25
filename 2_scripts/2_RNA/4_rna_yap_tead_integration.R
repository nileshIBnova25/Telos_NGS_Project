# WI-38 replicative senescence: RNA trend + YAP1/TEAD1 analysis
# Run from the Telos_NGS_Project repository root.

suppressPackageStartupMessages(library(readxl))

xlsx <- "1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx"
tpm_file <- "1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv"
outdir <- "1_data_raw/4_rnaseq/outputs/yap_tead_rna"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(xlsx)) stop("Workbook not found: ", xlsx)
if (!file.exists(tpm_file)) stop("TPM ground-truth table not found: ", tpm_file)

# 1. Read the existing DESeq2 time-course results.
de0 <- as.data.frame(read_excel(xlsx, sheet = "RS_deseq2"))
cat("RS_deseq2 columns:\n")
print(names(de0))

norm <- tolower(gsub("[^[:alnum:]]", "", names(de0)))

find_one <- function(index, label) {
    if (length(index) != 1L) {
        stop(
            "Could not identify exactly one ", label,
            ". Please share these column names: ",
            paste(names(de0), collapse = " | ")
        )
    }
    names(de0)[index]
}

gene_i <- which(norm %in% c(
    "gene", "genesymbol", "genename", "geneid", "geneidsymbol",
    "geneidorsymbol"
))
if (length(gene_i) != 1L) {
    gene_i <- grep("^gene($|symbol|name|id)", norm)
}
gene_col <- find_one(gene_i, "gene identifier column")

lfc_i <- which(norm %in% c(
    "log2foldchange", "log2fc", "logfc", "log2foldchangeestimate"
))
if (length(lfc_i) != 1L) {
    lfc_i <- grep("log2foldchange|^log2fc$", norm)
}
lfc_col <- find_one(lfc_i, "log2 fold-change column")

padj_i <- which(norm %in% c(
    "padj", "adjustedpvalue", "fdr", "fdradjustedpvalue",
    "padjusted", "padjvalue", "pvalueadjusted"
))
if (length(padj_i) != 1L) {
    padj_i <- grep("padj|adjustedpvalue|^fdr$", norm)
}
padj_col <- find_one(padj_i, "adjusted p-value column")

de <- data.frame(
    gene_symbol = trimws(as.character(de0[[gene_col]])),
    log2FC_time = suppressWarnings(as.numeric(as.character(de0[[lfc_col]]))),
    padj = suppressWarnings(as.numeric(as.character(de0[[padj_col]]))),
    stringsAsFactors = FALSE
)

de <- de[!is.na(de$gene_symbol) & nzchar(de$gene_symbol), ]
de$gene_key <- toupper(de$gene_symbol)

# Keep one entry per gene symbol, preferring the smallest adjusted p-value.
de <- de[order(is.na(de$padj), de$padj, -abs(de$log2FC_time), na.last = TRUE), ]
de <- de[!duplicated(de$gene_key), ]

# 2. Join the descriptive young-versus-late mean-TPM comparison.
tpm <- read.csv(tpm_file, stringsAsFactors = FALSE, check.names = FALSE)
required_tpm <- c(
    "gene_id_or_symbol", "young_tpm_mean", "rs_tpm_mean",
    "log2FC_rs_vs_young"
)
missing_tpm <- setdiff(required_tpm, names(tpm))
if (length(missing_tpm)) {
    stop("Missing columns in TPM table: ", paste(missing_tpm, collapse = ", "))
}

tpm$gene_key <- toupper(trimws(as.character(tpm$gene_id_or_symbol)))
tpm <- tpm[!duplicated(tpm$gene_key), ]
idx <- match(de$gene_key, tpm$gene_key)

de$young_tpm_mean <- tpm$young_tpm_mean[idx]
de$rs_tpm_mean <- tpm$rs_tpm_mean[idx]
de$log2FC_rs_vs_young <- tpm$log2FC_rs_vs_young[idx]

de$time_direction <- ifelse(
    is.na(de$log2FC_time), NA,
    ifelse(de$log2FC_time > 0, "increasing_over_time",
           ifelse(de$log2FC_time < 0, "decreasing_over_time", "zero"))
)

de$endpoint_direction <- ifelse(
    is.na(de$log2FC_rs_vs_young), NA,
    ifelse(de$log2FC_rs_vs_young > 0, "higher_mean_TPM_in_late",
           ifelse(de$log2FC_rs_vs_young < 0, "lower_mean_TPM_in_late", "zero"))
)

de$direction_concordance <- ifelse(
    is.na(de$log2FC_time) | is.na(de$log2FC_rs_vs_young),
    NA,
    ifelse(
        sign(de$log2FC_time) == sign(de$log2FC_rs_vs_young),
        "same_direction", "different_direction"
    )
)

write.csv(
    de,
    file.path(outdir, "RS_time_trends_with_endpoint_TPM.csv"),
    row.names = FALSE
)

# 3. Make statistical up/down lists from the provided time-course results.
# FDR < 0.01 is the paper's DE threshold.
# The stricter +/-0.5 log2FC lists mirror the gene-selection cutoffs
# used for the paper's Figure 6 regulator/motif analyses.
sig_up <- de[!is.na(de$padj) & de$padj < 0.01 &
             !is.na(de$log2FC_time) & de$log2FC_time > 0, ]
sig_down <- de[!is.na(de$padj) & de$padj < 0.01 &
               !is.na(de$log2FC_time) & de$log2FC_time < 0, ]

strict_up <- sig_up[sig_up$log2FC_time > 0.5, ]
strict_down <- sig_down[sig_down$log2FC_time < -0.5, ]

write.csv(sig_up, file.path(outdir, "RS_time_up_FDR01.csv"), row.names = FALSE)
write.csv(sig_down, file.path(outdir, "RS_time_down_FDR01.csv"), row.names = FALSE)
write.csv(strict_up, file.path(outdir, "RS_time_up_FDR01_log2FC_gt0.5.csv"), row.names = FALSE)
write.csv(strict_down, file.path(outdir, "RS_time_down_FDR01_log2FC_ltminus0.5.csv"), row.names = FALSE)

write_gene_list <- function(x, path) {
    genes <- unique(as.character(x$gene_symbol))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    writeLines(genes, path, useBytes = TRUE)
}

# Lists for g:Profiler/PANTHER (FDR < 0.01, split by direction).
write_gene_list(sig_up, file.path(outdir, "RS_time_up_FDR01.txt"))
write_gene_list(sig_down, file.path(outdir, "RS_time_down_FDR01.txt"))

# Lists for LISA (more restrictive effect-size cutoff; max 500 genes).
write_gene_list(head(strict_up, 500), file.path(outdir, "LISA_up_top500.txt"))
write_gene_list(head(strict_down, 500), file.path(outdir, "LISA_down_top500.txt"))

# Use all genes tested in RS_deseq2 as the enrichment background.
writeLines(
    sort(unique(de$gene_symbol)),
    file.path(outdir, "RS_all_tested_gene_background.txt"),
    useBytes = TRUE
)

# 4. Focused candidate panel based on genes discussed in the paper.
# These categories distinguish regulators, selected YAP/TEAD targets,
# and phenotype markers; they are not all direct targets.
panel <- data.frame(
    gene_symbol = c(
        "YAP1", "TEAD1", "TEAD4", "WWTR1",
        "CTGF", "CYR61", "TGFB2", "THBS1",
        "TOP2A", "CDC20", "BIRC5", "CDK9",
        "ACTA2", "FSTL1", "FBN1", "COL1A1", "COL1A2", "COL3A1"
    ),
    evidence_category = c(
        rep("Regulator_gene_expression_not_activity", 4),
        rep("YAP_TEAD_target_or_candidate_discussed_in_paper", 8),
        rep("Myofibroblast_or_collagen_marker_not_proof_of_direct_target", 6)
    ),
    stringsAsFactors = FALSE
)

pidx <- match(toupper(panel$gene_symbol), de$gene_key)
for (col in c(
    "log2FC_time", "padj", "young_tpm_mean", "rs_tpm_mean",
    "log2FC_rs_vs_young", "time_direction", "endpoint_direction",
    "direction_concordance"
)) {
    panel[[col]] <- de[[col]][pidx]
}

write.csv(
    panel,
    file.path(outdir, "YAP_TEAD_candidate_panel.csv"),
    row.names = FALSE
)

# 5. GSEA on the full time-course ranking using MSigDB Hallmark sets.
# This follows the paper's use of a full DESeq2 log2FC ranking + fgsea.
if (requireNamespace("fgsea", quietly = TRUE) &&
    requireNamespace("msigdbr", quietly = TRUE)) {

    msig <- tryCatch(
        msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
        error = function(e) {
            msigdbr::msigdbr(species = "Homo sapiens", category = "H")
        }
    )

    if (!all(c("gs_name", "gene_symbol") %in% names(msig))) {
        stop("Unexpected msigdbr columns: ", paste(names(msig), collapse = ", "))
    }

    pathways <- split(msig$gene_symbol, msig$gs_name)

    rank_df <- de[
        is.finite(de$log2FC_time) & !is.na(de$gene_symbol) &
        nzchar(de$gene_symbol), c("gene_symbol", "log2FC_time")
    ]
    rank_df <- rank_df[!duplicated(rank_df$gene_symbol), ]

    ranks <- setNames(rank_df$log2FC_time, rank_df$gene_symbol)
    ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))]
    ranks <- sort(ranks, decreasing = TRUE)

    set.seed(42)
    fg <- fgsea::fgsea(
        pathways = pathways,
        stats = ranks,
        minSize = 15,
        maxSize = 500,
        eps = 0
    )

    fg <- as.data.frame(fg)
    if ("leadingEdge" %in% names(fg)) {
        fg$leadingEdge <- vapply(
            fg$leadingEdge, function(x) paste(x, collapse = ";"),
            character(1)
        )
    }
    fg <- fg[order(fg$padj, na.last = TRUE), ]

    write.csv(
        fg,
        file.path(outdir, "GSEA_Hallmark_RS_time_trend.csv"),
        row.names = FALSE
    )

    # Plot the 20 most significant Hallmark sets, showing enrichment direction.
    top <- fg[!is.na(fg$padj), ]
    top <- head(top[order(top$padj), ], 20)

    if (nrow(top) > 0) {
        top <- top[order(top$NES), ]
        pdf(file.path(outdir, "GSEA_Hallmark_top20.pdf"), width = 10, height = 8)
        par(mar = c(5, 16, 3, 1))
        barplot(
            top$NES,
            names.arg = top$pathway,
            horiz = TRUE,
            las = 1,
            cex.names = 0.65,
            xlab = "Normalized enrichment score (NES)",
            main = "Hallmark GSEA: replicative-senescence time trend"
        )
        abline(v = 0, lty = 2)
        dev.off()
    }
} else {
    message(
        "GSEA skipped: install 'msigdbr' and Bioconductor 'fgsea', then rerun."
    )
}

cat("\nRNA/YAP-TEAD analysis completed.\n")
cat("RS_deseq2 genes read:", nrow(de), "\n")
cat("FDR < 0.01, increasing over time:", nrow(sig_up), "\n")
cat("FDR < 0.01, decreasing over time:", nrow(sig_down), "\n")
cat("Strict up list (log2FC > 0.5):", nrow(strict_up), "\n")
cat("Strict down list (log2FC < -0.5):", nrow(strict_down), "\n")
cat("Output directory:", outdir, "\n")
cat("Reminder: time trends are not automatically a direct PDL20-vs-late contrast.\n")
