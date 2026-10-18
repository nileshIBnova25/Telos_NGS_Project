from pathlib import Path

import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]

TRUTH_FILE = (
    ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)
PRED_FILE = (
    ROOT
    / "1_data_raw/3_alphagenome/02_atac_run/"
      "benchmark_10000_alphagenome_atac_pred.csv"
)
OUT_FILE = ROOT / "3_results/tables/prediction_senescence_correlations.csv"

measurements = [
    "young_mean",
    "rs_mean",
    "log2FC_rs_vs_young",
]

truth = pd.read_csv(TRUTH_FILE)
pred = pd.read_csv(PRED_FILE)

df = truth.merge(
    pred[["chr", "start", "end", "pred_atac_mean"]],
    on=["chr", "start", "end"],
    how="inner",
    validate="one_to_one",
)

df["pred_atac_mean"] = pd.to_numeric(
    df["pred_atac_mean"], errors="coerce"
)

classes = ["ALL", "RS_gained", "RS_lost", "Stable", "Negative"]
results = []

for label in classes:
    if label == "ALL":
        group = df.copy()
    else:
        group = df[df["label"] == label].copy()

    for measurement in measurements:
        valid = group[
            ["pred_atac_mean", measurement]
        ].dropna()

        n = len(valid)

        if n < 3:
            pearson_r, pearson_p = float("nan"), float("nan")
            spearman_r, spearman_p = float("nan"), float("nan")
        elif (
            valid["pred_atac_mean"].nunique() < 2
            or valid[measurement].nunique() < 2
        ):
            pearson_r, pearson_p = float("nan"), float("nan")
            spearman_r, spearman_p = float("nan"), float("nan")
        else:
            pearson_r, pearson_p = pearsonr(
                valid["pred_atac_mean"], valid[measurement]
            )
            spearman_r, spearman_p = spearmanr(
                valid["pred_atac_mean"], valid[measurement]
            )

        results.append({
            "class": label,
            "measurement": measurement,
            "n": n,
            "pearson_r": pearson_r,
            "pearson_p": pearson_p,
            "spearman_r": spearman_r,
            "spearman_p": spearman_p,
        })

result_df = pd.DataFrame(results)
OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
result_df.to_csv(OUT_FILE, index=False)

print("PREDICTION VS SENESCENCE MEASUREMENTS")
print("------------------------------------")
print(result_df.to_string(index=False, float_format=lambda x: f"{x:.6g}"))
print(f"\nSaved results to: {OUT_FILE}")
