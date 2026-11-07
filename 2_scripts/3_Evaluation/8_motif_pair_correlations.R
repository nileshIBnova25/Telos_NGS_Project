library(data.table)

root <- path.expand("~/Telos_NGS_Project")

input_file <- file.path(
  root, "3_results/tables/unified_motif_analysis",
  "all_nine_motifs_combined.csv"
)

ids_file <- file.path(
  root, "3_results/tables/unified_motif_analysis",
  "complete_nine_motif_region_ids.txt"
)

output_dir <- file.path(
  root, "3_results/tables/motif_cooccurrence_score_correlations"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read the existing calibrated motif scores
dat <- fread(input_file)
ids <- scan(ids_file, what = character(), quiet = TRUE)

required <- c("region_id", "label", "motif", "score_excess")
if (!all(required %in% names(dat))) {
  stop("Required columns are missing from the input CSV.")
}

dat <- dat[region_id %in% ids, ..required]

motifs <- sort(unique(dat$motif))
if (length(motifs) != 9) {
  stop("Expected 9 motifs, found ", length(motifs))
}

# Verify one score per region and motif
if (anyDuplicated(dat, by = c("region_id", "motif"))) {
  stop("Duplicate region_id/motif combinations detected.")
}

wide <- dcast(
  dat,
  region_id + label ~ motif,
  value.var = "score_excess"
)

if (nrow(wide) != length(ids)) {
  stop("The number of wide-table regions does not match the ID file.")
}

if (!all(motifs %in% names(wide))) {
  stop("One or more motif columns are missing.")
}

cat("Regions:", nrow(wide), "\n")
cat("Motifs:", length(motifs), "\n")
cat("Classes:\n")
print(wide[, .N, by = label])

# Calculate Spearman correlations for each unique pair
classes <- c("ALL", "RS_gained", "RS_lost", "Stable", "Negative")
pairs <- combn(motifs, 2, simplify = FALSE)
results <- list()
k <- 1

for (cl in classes) {
  subset <- if (cl == "ALL") wide else wide[label == cl]

  for (pair in pairs) {
    x <- subset[[pair[1]]]
    y <- subset[[pair[2]]]
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]

    test <- tryCatch(
      suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)),
      error = function(e) NULL
    )

    results[[k]] <- data.table(
      class = cl,
      motif1 = pair[1],
      motif2 = pair[2],
      n = length(x),
      spearman_rho = if (is.null(test)) NA_real_ else unname(test$estimate),
      p_value = if (is.null(test)) NA_real_ else test$p.value
    )
    k <- k + 1
  }
}

res <- rbindlist(results)
res[, BH_q_value := p.adjust(p_value, method = "BH"), by = class]
setorder(res, class, BH_q_value)

out_file <- file.path(output_dir, "motif_pair_spearman_correlations.csv")
fwrite(res, out_file)

cat("\nSaved results to:\n", out_file, "\n")
cat("Total pairwise tests:", nrow(res), "\n")
cat("\nTop 10 correlations by class (smallest BH-adjusted q-values):\n")

for (cl in classes) {
  cat("\n---", cl, "---\n")
  print(
    res[class == cl][order(BH_q_value)][1:min(10, .N)],
    digits = 4
  )
}
