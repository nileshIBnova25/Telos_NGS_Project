from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

root = Path.cwd()

pred_path = root / "1_data_raw/3_alphagenome/02_atac_run/benchmark_10000_alphagenome_atac_pred.csv"
truth_path = root / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
motif_path = root / "1_data_raw/5_motif_analysis/tead1_dinucleotide_calibration_regions.csv"

out_dir = root / "3_results/tables/tead1_alphagenome_error"
out_dir.mkdir(parents=True, exist_ok=True)

pred = pd.read_csv(pred_path)
truth = pd.read_csv(truth_path)
motif = pd.read_csv(motif_path)

# Construct a calibrated score: observed maximum minus shuffled-null mean.
motif["tead1_score_excess"] = (
    motif["observed_max_score"] - motif["null_mean_max_score"]
)

df = pred.merge(
    truth[["name", "rs_mean", "young_mean", "log2FC_rs_vs_young"]],
    left_on="name", right_on="name", how="inner",
    validate="one_to_one"
)

df = df.merge(
    motif[["region_id", "tead1_score_excess", "n_shuffles",
           "empirical_p", "BH_q_value"]],
    left_on="name", right_on="region_id", how="inner",
    validate="one_to_one"
)

df["signed_residual"] = df["pred_atac_mean"] - df["rs_mean"]
df["absolute_error"] = df["signed_residual"].abs()

df.to_csv(out_dir / "tead1_region_level_prediction_errors.csv", index=False)

print("Matched regions:", len(df))
print("Classes:")
print(df["label"].value_counts().to_string())
print("\nShuffles per region:")
print(df["n_shuffles"].value_counts().sort_index().to_string())

tests = []
groups = [("ALL", df)] + [
    (label, group) for label, group in df.groupby("label")
]

for group_name, group in groups:
    for outcome in ["absolute_error", "signed_residual"]:
        valid = group[["tead1_score_excess", outcome]].replace(
            [np.inf, -np.inf], np.nan
        ).dropna()

        rho, p = spearmanr(
            valid["tead1_score_excess"], valid[outcome]
        )

        tests.append({
            "group": group_name,
            "outcome": outcome,
            "n": len(valid),
            "spearman_rho": rho,
            "p_value": p
        })

results = pd.DataFrame(tests)
results["BH_q_value"] = multipletests(
    results["p_value"], method="fdr_bh"
)[1]

results.to_csv(
    out_dir / "tead1_error_correlations.csv", index=False
)

print("\nTEAD1 score excess versus AlphaGenome prediction error")
print(results.round(6).to_string(index=False))

class_summary = df.groupby("label").agg(
    n=("name", "size"),
    mean_tead1_score_excess=("tead1_score_excess", "mean"),
    median_tead1_score_excess=("tead1_score_excess", "median"),
    mean_absolute_error=("absolute_error", "mean"),
    median_absolute_error=("absolute_error", "median"),
    mean_signed_residual=("signed_residual", "mean")
).reset_index()

class_summary.to_csv(
    out_dir / "tead1_error_class_summary.csv", index=False
)

print("\nClass-level summary")
print(class_summary.round(5).to_string(index=False))

print("\nSaved results in:", out_dir)
