from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

root = Path.cwd()

truth = pd.read_csv(
    root / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)
cal = pd.read_csv(
    root / "1_data_raw/5_motif_analysis/tead1_dinucleotide_calibration_regions.csv"
)

cal["tead1_score_excess"] = (
    cal["observed_max_score"] - cal["null_mean_max_score"]
)

df = truth.merge(
    cal[["region_id", "tead1_score_excess", "n_shuffles"]],
    left_on="name",
    right_on="region_id",
    how="inner",
    validate="one_to_one"
)

df = df.dropna(
    subset=["tead1_score_excess", "log2FC_rs_vs_young"]
).copy()

print("Regions analysed:", len(df))
print("Classes:")
print(df["label"].value_counts().to_string())

rows = []

groups = [("ALL", df)] + list(df.groupby("label"))

for label, group in groups:
    rho, p = spearmanr(
        group["tead1_score_excess"],
        group["log2FC_rs_vs_young"]
    )

    rows.append({
        "label": label,
        "n": len(group),
        "spearman_rho": rho,
        "p_value": p
    })

results = pd.DataFrame(rows)
results["BH_q_value"] = multipletests(
    results["p_value"], method="fdr_bh"
)[1]

out_dir = root / "3_results/tables/tead1_accessibility_change"
out_dir.mkdir(parents=True, exist_ok=True)

results.to_csv(
    out_dir / "tead1_calibrated_score_accessibility_correlations.csv",
    index=False
)

print("\nCalibrated TEAD1 score excess versus observed ATAC-seq log2FC")
print(results.round(6).to_string(index=False))

print("\nSaved to:", out_dir)
