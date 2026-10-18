#!/usr/bin/env python3

"""Evaluate AlphaGenome ATAC signal separation across benchmark classes."""

from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd

from scipy.stats import kruskal, mannwhitneyu
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix,
)

ROOT = Path(__file__).resolve().parents[2]

GROUNDTRUTH = (
    ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)
PREDICTIONS = (
    ROOT
    / "1_data_raw/3_alphagenome/02_atac_run/"
    "benchmark_10000_alphagenome_atac_pred.csv"
)

OUT = ROOT / "3_results"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"

CLASSES = ["RS_gained", "RS_lost", "Stable", "Negative"]


def bh_adjust(p_values):
    """Benjamini-Hochberg false-discovery-rate adjustment."""
    p = np.asarray(p_values, dtype=float)
    adjusted = np.full(len(p), np.nan)

    valid_indices = np.where(np.isfinite(p))[0]
    if len(valid_indices) == 0:
        return adjusted

    valid_p = p[valid_indices]
    order = np.argsort(valid_p)
    sorted_p = valid_p[order]
    m = len(sorted_p)

    sorted_adjusted = sorted_p * m / np.arange(1, m + 1)
    sorted_adjusted = np.minimum.accumulate(sorted_adjusted[::-1])[::-1]
    sorted_adjusted = np.clip(sorted_adjusted, 0, 1)

    original_order = np.empty(m, dtype=int)
    original_order[order] = np.arange(m)
    adjusted[valid_indices] = sorted_adjusted[original_order]

    return adjusted


def main():
    if not GROUNDTRUTH.exists():
        raise FileNotFoundError(f"Missing file: {GROUNDTRUTH}")
    if not PREDICTIONS.exists():
        raise FileNotFoundError(f"Missing file: {PREDICTIONS}")

    truth = pd.read_csv(GROUNDTRUTH)
    predictions = pd.read_csv(PREDICTIONS)

    keys = ["chr", "start", "end"]
    required_truth = set(keys + ["label", "young_mean", "rs_mean"])
    required_pred = set(keys + ["pred_atac_mean"])

    if not required_truth.issubset(truth.columns):
        raise ValueError(
            f"Ground-truth columns missing: "
            f"{sorted(required_truth - set(truth.columns))}"
        )
    if not required_pred.issubset(predictions.columns):
        raise ValueError(
            f"Prediction columns missing: "
            f"{sorted(required_pred - set(predictions.columns))}"
        )

    truth["chr"] = truth["chr"].astype(str)
    predictions["chr"] = predictions["chr"].astype(str)

    for frame in (truth, predictions):
        frame["start"] = pd.to_numeric(frame["start"], errors="coerce")
        frame["end"] = pd.to_numeric(frame["end"], errors="coerce")

    if truth.duplicated(keys).any() or predictions.duplicated(keys).any():
        raise ValueError("Duplicate genomic coordinates found.")

    data = truth.merge(
        predictions[keys + ["pred_atac_mean"]],
        on=keys,
        how="inner",
        validate="one_to_one",
    )

    if len(data) != len(truth) or len(data) != len(predictions):
        raise ValueError(
            f"Coordinate mismatch: truth={len(truth)}, "
            f"predictions={len(predictions)}, matched={len(data)}"
        )

    data["pred_atac_mean"] = pd.to_numeric(
        data["pred_atac_mean"], errors="coerce"
    )
    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna(subset=["label", "pred_atac_mean"])

    observed_classes = set(data["label"].unique())
    if observed_classes != set(CLASSES):
        raise ValueError(
            f"Unexpected benchmark classes: {sorted(observed_classes)}"
        )

    # 1. Descriptive statistics for each benchmark class.
    summary = (
        data.groupby("label")["pred_atac_mean"]
        .agg(["count", "mean", "std", "median", "min", "max"])
        .reindex(CLASSES)
        .reset_index()
    )
    summary.to_csv(TABLES / "class_separation_summary.csv", index=False)

    # 2. Overall Kruskal-Wallis test.
    groups = [
        data.loc[data["label"] == label, "pred_atac_mean"].to_numpy()
        for label in CLASSES
    ]
    h_stat, kw_p = kruskal(*groups)

    pd.DataFrame([{
        "test": "Kruskal-Wallis",
        "classes": ";".join(CLASSES),
        "n": len(data),
        "statistic_H": h_stat,
        "p_value": kw_p,
    }]).to_csv(TABLES / "class_separation_kruskal_wallis.csv", index=False)

    # 3. Pairwise Mann-Whitney U tests with BH correction.
    pair_rows = []
    for class_a, class_b in combinations(CLASSES, 2):
        a = data.loc[data["label"] == class_a, "pred_atac_mean"]
        b = data.loc[data["label"] == class_b, "pred_atac_mean"]

        test = mannwhitneyu(a, b, alternative="two-sided")
        pair_rows.append({
            "class_a": class_a,
            "class_b": class_b,
            "n_a": len(a),
            "n_b": len(b),
            "median_a": a.median(),
            "median_b": b.median(),
            "median_difference_a_minus_b": a.median() - b.median(),
            "U_statistic": test.statistic,
            "p_value": test.pvalue,
        })

    pairwise = pd.DataFrame(pair_rows)
    pairwise["p_adjusted_BH"] = bh_adjust(pairwise["p_value"])
    pairwise.to_csv(TABLES / "class_separation_pairwise.csv", index=False)

    # 4. One-vs-rest exploratory classification using predicted signal
    # as the score. These metrics do not imply passage-aware prediction.
    classification_rows = []
    for label in CLASSES:
        y_true = (data["label"] == label).astype(int).to_numpy()
        scores = data["pred_atac_mean"].to_numpy()

        # A constant score cannot rank positives above negatives.
        auc = roc_auc_score(y_true, scores)
        ap = average_precision_score(y_true, scores)

        # Illustrative threshold only; not optimized on this dataset.
        threshold = data["pred_atac_mean"].median()
        y_pred = (scores >= threshold).astype(int)

        tn, fp, fn, tp = confusion_matrix(
            y_true, y_pred, labels=[0, 1]
        ).ravel()

        classification_rows.append({
            "target_class": label,
            "n_positive": int(y_true.sum()),
            "n_negative": int((1 - y_true).sum()),
            "roc_auc": auc,
            "average_precision": ap,
            "threshold_median_all_predictions": threshold,
            "precision_at_threshold": precision_score(
                y_true, y_pred, zero_division=0
            ),
            "recall_at_threshold": recall_score(
                y_true, y_pred, zero_division=0
            ),
            "f1_at_threshold": f1_score(
                y_true, y_pred, zero_division=0
            ),
            "true_positive": int(tp),
            "false_positive": int(fp),
            "true_negative": int(tn),
            "false_negative": int(fn),
        })

    pd.DataFrame(classification_rows).to_csv(
        TABLES / "class_separation_one_vs_rest.csv", index=False
    )

    print("Class-separation evaluation completed.")
    print(f"Matched regions evaluated: {len(data)}")
    print(f"Results directory: {TABLES}")
    print("\nClass summary:")
    print(summary.to_string(index=False))
    print("\nKruskal-Wallis test:")
    print(f"H = {h_stat:.4f}, p-value = {kw_p:.6g}")
    print("\nPairwise comparisons (BH-adjusted p-values):")
    print(
        pairwise[
            ["class_a", "class_b", "p_value", "p_adjusted_BH"]
        ].to_string(index=False)
    )
    print("\nOne-vs-rest ROC-AUC and average precision:")
    print(
        pd.DataFrame(classification_rows)[
            ["target_class", "roc_auc", "average_precision"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    main()
