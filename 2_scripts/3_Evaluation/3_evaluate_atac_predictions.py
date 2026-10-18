#!/usr/bin/env python3

"""Evaluate AlphaGenome ATAC predictions against empirical benchmark signals."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error


ROOT = Path(__file__).resolve().parents[2]

GROUNDTRUTH = (
    ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)

PREDICTIONS = (
    ROOT
    / "1_data_raw/3_alphagenome/02_atac_run/"
    "benchmark_10000_alphagenome_atac_pred.csv"
)

OUTPUT_DIR = ROOT / "3_results/tables"
OUTPUT_FILE = OUTPUT_DIR / "alphagenome_atac_evaluation.csv"


def correlation(x, y, method):
    """Calculate a correlation, returning NaN if it cannot be estimated."""
    if len(x) < 2 or x.nunique() < 2 or y.nunique() < 2:
        return np.nan

    if method == "pearson":
        return pearsonr(x, y).statistic

    return spearmanr(x, y).statistic


def main():
    if not GROUNDTRUTH.exists():
        raise FileNotFoundError(f"Ground truth not found: {GROUNDTRUTH}")

    if not PREDICTIONS.exists():
        raise FileNotFoundError(f"Prediction file not found: {PREDICTIONS}")

    truth = pd.read_csv(GROUNDTRUTH)
    predictions = pd.read_csv(PREDICTIONS)

    required_truth = {
        "chr", "start", "end", "name", "label",
        "young_mean", "rs_mean", "log2FC_rs_vs_young",
    }
    required_predictions = {
        "chr", "start", "end", "pred_atac_mean",
    }

    missing_truth = required_truth - set(truth.columns)
    missing_predictions = required_predictions - set(predictions.columns)

    if missing_truth:
        raise ValueError(f"Missing ground-truth columns: {sorted(missing_truth)}")

    if missing_predictions:
        raise ValueError(
            f"Missing prediction columns: {sorted(missing_predictions)}"
        )

    # Match predictions to ground truth by genomic coordinates.
    keys = ["chr", "start", "end"]

    truth["chr"] = truth["chr"].astype(str)
    predictions["chr"] = predictions["chr"].astype(str)

    for frame in (truth, predictions):
        frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
        frame["end"] = pd.to_numeric(frame["end"], errors="coerce")

    if truth.duplicated(keys).any():
        raise ValueError("Ground truth contains duplicate genomic coordinates.")

    if predictions.duplicated(keys).any():
        raise ValueError("Predictions contain duplicate genomic coordinates.")

    merged = truth.merge(
        predictions[keys + ["pred_atac_mean"]],
        on=keys,
        how="inner",
        validate="one_to_one",
    )

    if len(merged) != len(truth) or len(merged) != len(predictions):
        raise ValueError(
            "Ground truth and predictions do not match one-to-one: "
            f"ground truth={len(truth)}, predictions={len(predictions)}, "
            f"matched={len(merged)}."
        )

    numeric_columns = [
        "young_mean",
        "rs_mean",
        "log2FC_rs_vs_young",
        "pred_atac_mean",
    ]

    for column in numeric_columns:
        merged[column] = pd.to_numeric(merged[column], errors="coerce")

    # Overall prediction-to-signal agreement.
    # AlphaGenome's prediction is not a separate prediction for each cell state.
    comparisons = [
        ("pred_vs_young", "pred_atac_mean", "young_mean"),
        ("pred_vs_rs", "pred_atac_mean", "rs_mean"),
    ]

    rows = []

    for comparison_name, prediction_col, target_col in comparisons:
        valid = merged[[prediction_col, target_col]].replace(
            [np.inf, -np.inf], np.nan
        ).dropna()

        x = valid[prediction_col]
        y = valid[target_col]

        rows.append({
            "analysis": comparison_name,
            "label": "ALL",
            "n": len(valid),
            "pearson_r": correlation(x, y, "pearson"),
            "spearman_r": correlation(x, y, "spearman"),
            "mae": mean_absolute_error(y, x) if len(valid) else np.nan,
            "rmse": (
                np.sqrt(mean_squared_error(y, x))
                if len(valid) else np.nan
            ),
        })

    # Repeat the comparisons within each benchmark class.
    for label, group in merged.groupby("label", sort=True):
        for comparison_name, prediction_col, target_col in comparisons:
            valid = group[[prediction_col, target_col]].replace(
                [np.inf, -np.inf], np.nan
            ).dropna()

            x = valid[prediction_col]
            y = valid[target_col]

            rows.append({
                "analysis": comparison_name,
                "label": label,
                "n": len(valid),
                "pearson_r": correlation(x, y, "pearson"),
                "spearman_r": correlation(x, y, "spearman"),
                "mae": mean_absolute_error(y, x) if len(valid) else np.nan,
                "rmse": (
                    np.sqrt(mean_squared_error(y, x))
                    if len(valid) else np.nan
                ),
            })

    # Descriptive benchmark-class summaries.
    class_summary = (
        merged.groupby("label", sort=True)
        .agg(
            n=("name", "size"),
            young_mean_avg=("young_mean", "mean"),
            rs_mean_avg=("rs_mean", "mean"),
            log2FC_avg=("log2FC_rs_vs_young", "mean"),
            predicted_atac_avg=("pred_atac_mean", "mean"),
        )
        .reset_index()
    )

    for _, row in class_summary.iterrows():
        rows.append({
            "analysis": "class_summary",
            "label": row["label"],
            "n": int(row["n"]),
            "pearson_r": np.nan,
            "spearman_r": np.nan,
            "mae": np.nan,
            "rmse": np.nan,
        })

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT_FILE, index=False)

    print("AlphaGenome ATAC evaluation script completed.")
    print(f"Ground-truth regions: {len(truth)}")
    print(f"Prediction rows:      {len(predictions)}")
    print(f"Matched regions:      {len(merged)}")
    print(f"Benchmark classes:    {merged['label'].nunique()}")
    print(f"Results saved to:     {OUTPUT_FILE}")
    print("\nCorrelation and error metrics:")
    print(
        results[results["analysis"].isin(["pred_vs_young", "pred_vs_rs"])]
        .to_string(index=False)
    )
    print("\nClass-level descriptive summary:")
    print(class_summary.to_string(index=False))


if __name__ == "__main__":
    main()
