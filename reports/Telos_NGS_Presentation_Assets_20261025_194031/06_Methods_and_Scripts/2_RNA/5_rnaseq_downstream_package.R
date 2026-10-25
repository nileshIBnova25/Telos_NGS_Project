#!/usr/bin/env Rscript

# WI-38 RNA-seq downstream package
# Run from the Telos_NGS_Project root:
# Rscript 2_scripts/2_RNA/5_rnaseq_downstream_package.R
#
# This script does not rerun DESeq2. It audits the study-provided RS_deseq2
# table, rebuilds clean gene lists/background, and creates descriptive figures.

options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c("readxl", "dplyr", "tidyr", "ggplot2", "readr", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste(sprintf("'%s'", missing_packages), collapse = ", "), "))")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

root <- getwd()
workbook <- file.path(root, "1_data_raw", "GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx")
existing_out <- file.path(root, "1_data_raw", "4_rnaseq", "outputs")
existing_trend <- file.path(existing_out, "yap_tead_rna", "RS_time_trends_with_endpoint_TPM.csv")
gsea_file <- file.path(existing_out, "yap_tead_rna", "GSEA_Hallmark_RS_time_trend.csv")
alpha_file <- file.path(existing_out, "alphagenome", "wi38_rnaseq_alphagenome_integrated.csv")
out_dir <- file.path(existing_out, "rnaseq_downstream_package")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(workbook)) stop("Workbook not found: ", workbook, "\nRun this script from the project root.")

# Helpers -------------------------------------------------------------------
first_matching_col <- function(nms, candidates) {
  hit <- candidates[tolower(candidates) %in% tolower(nms)]
  if (length(hit)) return(nms[match(tolower(hit[1]), tolower(nms))])
  NA_character_
}

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  # Prevent date strings from being silently treated as gene symbols.
  x <- x[!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)]
  unique(x)
}

write_gene_list <- function(x, path) {
  x <- clean_symbols(x)
  writeLines(x, path, useBytes = TRUE)
  invisible(length(x))
}

save_plot <- function(p, name, width = 10, height = 7) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = "pdf")
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = 300)
}

# 1. Read original study tables ------------------------------------------------
des <- readxl::read_excel(workbook, sheet = "RS_deseq2", guess_max = 100000) |> as.data.frame()
raw_tpm <- readxl::read_excel(workbook, sheet = "RS_TC_raw_TPM", guess_max = 100000) |> as.data.frame()

gene_col_des <- first_matching_col(names(des), c("gene", "gene_id_or_symbol", "symbol", "gene_symbol"))
lfc_col <- first_matching_col(names(des), c("log2FoldChange", "log2FC", "log2FC_time", "log2foldchange"))
padj_col <- first_matching_col(names(des), c("padj", "FDR", "adj.P.Val", "adjusted_p_value"))
pval_col <- first_matching_col(names(des), c("pvalue", "pvalue", "P.Value"))
if (any(is.na(c(gene_col_des, lfc_col, padj_col)))) {
  stop("Could not identify required columns in RS_deseq2. Columns found: ", paste(names(des), collapse = ", "))
}

# Preserve the original column names, but standardise what we need downstream.
des2 <- data.frame(
  gene = trimws(as.character(des[[gene_col_des]])),
  baseMean = if ("baseMean" %in% names(des)) suppressWarnings(as.numeric(des$baseMean)) else NA_real_,
  log2FC_time = suppressWarnings(as.numeric(des[[lfc_col]])),
  pvalue = if (!is.na(pval_col)) suppressWarnings(as.numeric(des[[pval_col]])) else NA_real_,
  padj = suppressWarnings(as.numeric(des[[padj_col]])),
  stringsAsFactors = FALSE
)
des2 <- des2[!is.na(des2$gene) & nzchar(des2$gene), , drop = FALSE]
des2 <- des2[!grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", des2$gene), , drop = FALSE]
# If the source contains duplicate symbols, retain the row with the smallest adjusted p-value.
des2 <- des2[order(ifelse(is.na(des2$padj), Inf, des2$padj), -abs(des2$log2FC_time)), , drop = FALSE]
des2 <- des2[!duplicated(des2$gene), , drop = FALSE]

if (nrow(des2) < 1000) stop("Only ", nrow(des2), " usable gene identifiers found in RS_deseq2 after removing dates. Stop: inspect the workbook before proceeding.")

# Build background from the original DESeq2 tested-gene table, not the previously
# generated background file (which was observed to start with date-like values).
background_genes <- clean_symbols(des2$gene)
date_like_background <- sum(grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", background_genes))
if (length(background_genes) < 1000 || date_like_background > 0) {
  stop("Background validation failed. Valid unique identifiers: ", length(background_genes),
       ". Date-like identifiers remaining: ", date_like_background,
       ". Do not use this background for enrichment until checked.")
}
write_gene_list(background_genes, file.path(out_dir, "RS_all_tested_gene_background_FIXED.txt"))

up <- des2 |> filter(!is.na(padj), padj < 0.01, log2FC_time > 0) |> arrange(padj, desc(log2FC_time))
down <- des2 |> filter(!is.na(padj), padj < 0.01, log2FC_time < 0) |> arrange(padj, log2FC_time)
up_strict <- up |> filter(log2FC_time > 0.5)
down_strict <- down |> filter(log2FC_time < -0.5)

write_csv(up, file.path(out_dir, "RS_time_up_FDR01_FIXED.csv"))
write_csv(down, file.path(out_dir, "RS_time_down_FDR01_FIXED.csv"))
write_csv(up_strict, file.path(out_dir, "RS_time_up_FDR01_log2FC_gt0.5_FIXED.csv"))
write_csv(down_strict, file.path(out_dir, "RS_time_down_FDR01_log2FC_ltminus0.5_FIXED.csv"))
write_gene_list(up$gene, file.path(out_dir, "RS_time_up_FDR01_FIXED.txt"))
write_gene_list(down$gene, file.path(out_dir, "RS_time_down_FDR01_FIXED.txt"))
# Lists for LISA are limited to 450 genes to remain safely below the commonly stated <500 limit.
write_gene_list(head(up$gene, 450), file.path(out_dir, "LISA_up_top450_FIXED.txt"))
write_gene_list(head(down$gene, 450), file.path(out_dir, "LISA_down_top450_FIXED.txt"))

# 2. Read endpoint TPM comparison ---------------------------------------------
endpoint <- NULL
if (file.exists(existing_trend)) {
  endpoint <- read.csv(existing_trend, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  gt_file <- file.path(existing_out, "rnaseq_groundtruth_wi38_young_vs_rs.csv")
  if (file.exists(gt_file)) endpoint <- read.csv(gt_file, check.names = FALSE, stringsAsFactors = FALSE)
}
if (!is.null(endpoint)) {
  gene_col_endpoint <- first_matching_col(names(endpoint), c("gene_id_or_symbol", "gene", "symbol", "gene_symbol"))
  if (is.na(gene_col_endpoint)) {
    warning("Could not identify gene column in endpoint TPM table. Endpoint comparisons will be skipped.")
    endpoint <- NULL
  } else {
    names(endpoint)[names(endpoint) == gene_col_endpoint] <- "gene"
    endpoint$gene <- trimws(as.character(endpoint$gene))
    if (!"log2FC_rs_vs_young" %in% names(endpoint)) {
      if (all(c("young_tpm_mean", "rs_tpm_mean") %in% names(endpoint))) {
        endpoint$log2FC_rs_vs_young <- log2((as.numeric(endpoint$rs_tpm_mean) + 0.001) /
                                             (as.numeric(endpoint$young_tpm_mean) + 0.001))
      }
    }
    if (all(c("young_tpm_mean", "rs_tpm_mean") %in% names(endpoint))) {
      endpoint$young_tpm_mean <- as.numeric(endpoint$young_tpm_mean)
      endpoint$rs_tpm_mean <- as.numeric(endpoint$rs_tpm_mean)
    }
  }
}

# 3. Candidate gene panel: paper-linked target candidates plus senescence markers.
candidate_panel <- data.frame(
  gene = c("YAP1", "TEAD1", "TEAD4", "WWTR1",
           "CTGF", "CYR61", "TGFB2", "THBS1",
           "TOP2A", "CDC20", "BIRC5", "CDK9",
           "ACTA2", "FSTL1", "FBN1", "COL1A1", "COL1A2", "COL3A1",
           "CDKN1A", "CDKN2A", "GLB1", "NNMT"),
  category = c(rep("YAP/TEAD regulators", 4),
               rep("YAP/TEAD-related candidates", 4),
               rep("Cell-cycle-associated candidates", 4),
               rep("ECM/myofibroblast markers", 6),
               rep("Senescence markers", 4)),
  stringsAsFactors = FALSE
)
panel <- candidate_panel |> left_join(des2, by = "gene")
if (!is.null(endpoint)) panel <- panel |> left_join(endpoint, by = "gene", suffix = c("", "_endpoint"))
panel <- panel |>
  mutate(
    trend_significance = case_when(is.na(padj) ~ "Not present / no result",
                                   padj < 0.01 & log2FC_time > 0 ~ "Increasing over time (FDR < 0.01)",
                                   padj < 0.01 & log2FC_time < 0 ~ "Decreasing over time (FDR < 0.01)",
                                   TRUE ~ "Not significant at FDR < 0.01"),
    endpoint_direction = if ("log2FC_rs_vs_young" %in% names(panel))
      case_when(is.na(log2FC_rs_vs_young) ~ "Endpoint not available",
                log2FC_rs_vs_young > 0.25 ~ "Higher in late TPM group",
                log2FC_rs_vs_young < -0.25 ~ "Lower in late TPM group",
                TRUE ~ "Near-stable endpoint TPM") else "Endpoint not available",
    direction_concordance = if ("log2FC_rs_vs_young" %in% names(panel))
      case_when(is.na(log2FC_rs_vs_young) | is.na(log2FC_time) ~ "Not assessable",
                log2FC_time * log2FC_rs_vs_young > 0 ~ "Same direction (descriptive)",
                log2FC_time * log2FC_rs_vs_young < 0 ~ "Opposite direction (descriptive)",
                TRUE ~ "One direction near zero") else "Not assessable"
  )
write_csv(panel, file.path(out_dir, "candidate_gene_panel_summary.csv"))

# 4. Candidate gene expression trajectories over all recorded PDLs ----------------
gene_col_tpm <- first_matching_col(names(raw_tpm), c("gene", "gene_id_or_symbol", "symbol", "gene_symbol"))
sample_cols <- grep("^PDL[0-9]+", names(raw_tpm), value = TRUE)
trajectory_data <- NULL
if (is.na(gene_col_tpm) || length(sample_cols) < 3) {
  warning("Could not identify gene/sample columns in RS_TC_raw_TPM; trajectory plot skipped.")
} else {
  tpm2 <- raw_tpm |> mutate(.gene = trimws(as.character(.data[[gene_col_tpm]]))) |>
    select(.gene, all_of(sample_cols))
  trajectory_data <- tpm2 |>
    pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "TPM") |>
    mutate(TPM = suppressWarnings(as.numeric(as.character(TPM))),
           PDL = suppressWarnings(as.integer(sub("^PDL([0-9]+).*", "\\1", sample))),
           gene = .gene) |>
    filter(gene %in% candidate_panel$gene, is.finite(PDL), !is.na(TPM)) |>
    select(gene, sample, PDL, TPM) |>
    left_join(candidate_panel, by = "gene")
  write_csv(trajectory_data, file.path(out_dir, "candidate_gene_trajectories.csv"))
  if (nrow(trajectory_data) > 0) {
    p_traj <- ggplot(trajectory_data, aes(x = PDL, y = log2(TPM + 0.1))) +
      geom_point(alpha = 0.55, size = 1.2, position = position_jitter(width = 0.25, height = 0)) +
      stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 0.7, color = "black") +
      facet_wrap(~ gene, scales = "free_y", ncol = 4) +
      scale_x_continuous(breaks = sort(unique(trajectory_data$PDL))) +
      labs(title = "Expression trajectories of selected senescence-related genes",
           subtitle = "Raw TPM across the recorded WI-38 replicative-senescence time course; y-axis is log2(TPM + 0.1)",
           x = "Population doubling level (PDL)", y = "log2(TPM + 0.1)") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "bold", size = 8),
            plot.title = element_text(face = "bold"))
    save_plot(p_traj, "figure_02_candidate_gene_trajectories", width = 13, height = 11)
  }
}

# 5. GSEA directional summary -------------------------------------------------
if (file.exists(gsea_file)) {
  gsea <- read.csv(gsea_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c("pathway", "padj", "NES") %in% names(gsea))) {
    gsea$padj <- as.numeric(gsea$padj)
    gsea$NES <- as.numeric(gsea$NES)
    sig <- gsea |> filter(!is.na(padj), padj < 0.05, is.finite(NES))
    chosen <- bind_rows(
      sig |> filter(NES > 0) |> arrange(padj) |> slice_head(n = 10),
      sig |> filter(NES < 0) |> arrange(padj) |> slice_head(n = 10)
    ) |> distinct(pathway, .keep_all = TRUE) |>
      mutate(direction = ifelse(NES > 0, "Enriched toward increasing genes", "Enriched toward decreasing genes"),
             term = gsub("_", " ", sub("^HALLMARK_", "", pathway)),
             term = factor(term, levels = term[order(NES)]))
    write_csv(chosen |> arrange(padj), file.path(out_dir, "GSEA_Hallmark_top_directional_terms.csv"))
    if (nrow(chosen) > 0) {
      p_gsea <- ggplot(chosen, aes(x = NES, y = term, fill = direction)) +
        geom_col(width = 0.72) +
        geom_vline(xintercept = 0, linewidth = 0.4) +
        labs(title = "Hallmark GSEA: strongest directional enrichments",
             subtitle = "Ranked by the supplied DESeq2 time-trend statistic; only adjusted p < 0.05 shown",
             x = "Normalized enrichment score (NES)", y = NULL, fill = NULL) +
        theme_bw(base_size = 10) +
        theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
      save_plot(p_gsea, "figure_01_GSEA_Hallmark_directional", width = 11, height = 7.5)
    }
  } else warning("GSEA file did not contain pathway, padj, and NES columns; plot skipped.")
} else warning("GSEA file not found; GSEA figure skipped.")

# 6. Time-trend coefficient vs endpoint TPM direction ------------------------
if (!is.null(endpoint) && all(c("gene", "log2FC_rs_vs_young") %in% names(endpoint))) {
  comparison <- des2 |> inner_join(endpoint |> select(gene, log2FC_rs_vs_young,
                            any_of(c("young_tpm_mean", "rs_tpm_mean"))), by = "gene")
  comparison$max_endpoint_tpm <- if (all(c("young_tpm_mean", "rs_tpm_mean") %in% names(comparison)))
    pmax(comparison$young_tpm_mean, comparison$rs_tpm_mean, na.rm = TRUE) else Inf
  # For the correlation/plot, focus on genes with at least 1 TPM in one endpoint.
  comparison_plot <- comparison |> filter(is.finite(log2FC_time), is.finite(log2FC_rs_vs_young), max_endpoint_tpm >= 1)
  if (nrow(comparison_plot) > 20) {
    rho <- suppressWarnings(cor(comparison_plot$log2FC_time, comparison_plot$log2FC_rs_vs_young,
                                method = "spearman", use = "complete.obs"))
    direction_agreement <- mean(sign(comparison_plot$log2FC_time) == sign(comparison_plot$log2FC_rs_vs_young), na.rm = TRUE)
    agreement_table <- data.frame(
      metric = c("Spearman correlation: DESeq2 time trend vs descriptive endpoint log2FC",
                 "Direction agreement fraction among genes with max endpoint TPM >= 1",
                 "Genes included in comparison"),
      value = c(rho, direction_agreement, nrow(comparison_plot))
    )
    write_csv(agreement_table, file.path(out_dir, "time_trend_vs_endpoint_agreement_metrics.csv"))
    p_comp <- ggplot(comparison_plot, aes(x = log2FC_rs_vs_young, y = log2FC_time)) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.6) +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.6) +
      geom_point(alpha = 0.22, size = 0.7) +
      annotate("text", x = Inf, y = Inf,
               label = sprintf("Spearman rho = %.2f\nDirection agreement = %.1f%%", rho, 100 * direction_agreement),
               hjust = 1.05, vjust = 1.2, size = 3.4) +
      coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
      labs(title = "Study time-trend effects versus descriptive endpoint expression changes",
           subtitle = "Gene filter: at least 1 TPM in young or late endpoint; axes truncated visually to ±4",
           x = "Descriptive endpoint log2FC (late TPM group vs young TPM group)",
           y = "Study-provided DESeq2 time-trend log2FoldChange") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"))
    save_plot(p_comp, "figure_03_endpoint_vs_time_trend", width = 8.5, height = 7)
    write_csv(comparison, file.path(out_dir, "all_genes_time_trend_vs_endpoint_comparison.csv"))
  } else warning("Too few genes available for endpoint/time-trend comparison.")
}

# 7. AlphaGenome predicted RNA track vs observed average TPM (exploratory only) --
# This tests cross-gene concordance of absolute signal/expression, not senescence-specific prediction.
alpha_status <- "AlphaGenome integrated CSV not found; RNA-signal comparison skipped."
if (file.exists(alpha_file)) {
  ag <- read.csv(alpha_file, stringsAsFactors = FALSE, check.names = FALSE)
  ag_gene_col <- first_matching_col(names(ag), c("gene_id_or_symbol", "gene", "symbol", "gene_symbol"))
  pred_candidates <- grep("pred.*rna_seq.*promoter.*mean", names(ag), ignore.case = TRUE, value = TRUE)
  observed_young <- first_matching_col(names(ag), c("young_tpm_mean"))
  observed_late <- first_matching_col(names(ag), c("rs_tpm_mean"))
  if (!is.na(ag_gene_col) && length(pred_candidates) > 0 && !is.na(observed_young) && !is.na(observed_late)) {
    pred_col <- pred_candidates[1]
    ag2 <- data.frame(
      gene = trimws(as.character(ag[[ag_gene_col]])),
      pred_rna_signal = suppressWarnings(as.numeric(ag[[pred_col]])),
      young_tpm_mean = suppressWarnings(as.numeric(ag[[observed_young]])),
      rs_tpm_mean = suppressWarnings(as.numeric(ag[[observed_late]])),
      stringsAsFactors = FALSE
    ) |>
      mutate(mean_observed_tpm = (young_tpm_mean + rs_tpm_mean) / 2,
             log2_mean_observed_tpm = log2(mean_observed_tpm + 0.1)) |>
      filter(!is.na(gene), nzchar(gene), is.finite(pred_rna_signal), is.finite(log2_mean_observed_tpm)) |>
      distinct(gene, .keep_all = TRUE)
    if (nrow(ag2) >= 20 && length(unique(ag2$pred_rna_signal)) > 2) {
      rho_ag <- suppressWarnings(cor(ag2$pred_rna_signal, ag2$log2_mean_observed_tpm,
                                     method = "spearman", use = "complete.obs"))
      alpha_metrics <- data.frame(
        predicted_column = pred_col,
        n_genes = nrow(ag2),
        spearman_rho_predicted_promoter_RNA_signal_vs_log2_mean_TPM = rho_ag,
        interpretation = "Exploratory cross-gene association only; the predicted track is not a predicted young-to-senescent fold change."
      )
      write_csv(alpha_metrics, file.path(out_dir, "alphagenome_predRNA_vs_observed_TPM_metrics.csv"))
      write_csv(ag2, file.path(out_dir, "alphagenome_predRNA_vs_observed_TPM_data.csv"))
      p_ag <- ggplot(ag2, aes(x = pred_rna_signal, y = log2_mean_observed_tpm)) +
        geom_point(alpha = 0.45, size = 1) +
        labs(title = "AlphaGenome promoter RNA signal vs observed average TPM",
             subtitle = paste0("Exploratory across genes; Spearman rho = ", round(rho_ag, 3),
                               ". Not a prediction of the senescence-associated fold change."),
             x = paste0("AlphaGenome predicted promoter RNA signal (", pred_col, ")"),
             y = "log2(mean TPM across young and late groups + 0.1)") +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(face = "bold"))
      save_plot(p_ag, "figure_04_alphagenome_predRNA_vs_observed_TPM", width = 8.5, height = 7)
      alpha_status <- paste0("AlphaGenome comparison created using column: ", pred_col,
                             ". It is a cross-gene signal/expression association, not condition-specific validation.")
    } else {
      alpha_status <- "AlphaGenome columns were found, but insufficient variable valid predictions were available for a correlation."
    }
  } else {
    alpha_status <- paste0("AlphaGenome integrated file was found, but required columns were not detected. Columns: ", paste(names(ag), collapse = ", "))
  }
}

# 8. Validation summary and LISA handoff template ------------------------------
summary_lines <- c(
  "WI-38 RNA-seq downstream package",
  paste0("Workbook: ", workbook),
  paste0("DESeq2 source sheet: RS_deseq2; gene column = ", gene_col_des,
         "; LFC column = ", lfc_col, "; adjusted p-value column = ", padj_col),
  paste0("Unique tested genes after cleaning: ", length(background_genes)),
  paste0("Increasing time trend, FDR < 0.01: ", nrow(up)),
  paste0("Decreasing time trend, FDR < 0.01: ", nrow(down)),
  paste0("Increasing time trend, FDR < 0.01 and log2FC > 0.5: ", nrow(up_strict)),
  paste0("Decreasing time trend, FDR < 0.01 and log2FC < -0.5: ", nrow(down_strict)),
  paste0("Background file: ", file.path(out_dir, "RS_all_tested_gene_background_FIXED.txt")),
  paste0("Raw TPM gene column = ", gene_col_tpm, "; detected PDL sample columns = ", length(sample_cols)),
  paste0("Candidate trajectory rows: ", if (is.null(trajectory_data)) 0 else nrow(trajectory_data)),
  paste0("AlphaGenome audit: ", alpha_status),
  "Caution: RS_deseq2 is a study-provided time-trend model. It is not automatically a direct PDL20 vs selected late TPM contrast.",
  "Caution: TPM fold changes in this package are descriptive and should not replace count-based differential-expression inference.",
  "Caution: AlphaGenome predicted RNA tracks are sequence-based predictions for a selected model track/biosample; they are not observed WI-38 young/late expression tracks."
)
writeLines(summary_lines, file.path(out_dir, "RUN_SUMMARY.txt"))

lisa_template <- data.frame(
  gene_set = c("Increasing_over_time", "Increasing_over_time", "Decreasing_over_time", "Decreasing_over_time"),
  regulator = c("YAP1", "TEAD1", "YAP1", "TEAD1"),
  rank = NA_integer_,
  p_value = NA_real_,
  adjusted_p_value = NA_real_,
  dataset_or_track_label_from_LISA = NA_character_,
  notes = "Fill from the downloaded LISA detailed results; do not infer rank/p-value from plot labels alone.",
  stringsAsFactors = FALSE
)
write_csv(lisa_template, file.path(out_dir, "LISA_YAP1_TEAD1_result_extraction_template.csv"))

cat(paste(summary_lines, collapse = "\n"), "\n\nOutput directory:\n", out_dir, "\n")
