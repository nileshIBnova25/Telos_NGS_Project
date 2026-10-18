#!/usr/bin/env python3

"""Follow-up evaluation: ranking direction, effect sizes, and baseline adjustment."""

from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score


ROOT = Path(__file__).resolve().parents[2]

GROUNDTRUTH = (
    ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)
PREDICTIONS = (
    ROOT
    / "1_data_raw/3_alphagenome/02_atac_run/"
    "benchmark_10000_alphagenome_atac_pred.csv"
)
OUT = ROOT / "3_results/tables"

CLASSES = ["RS_gained", "RS_lost", "Stable", "Negative"]


def cliff_delta(a, b):
    """Effect size: probability(a > b) minus probability(a < b)."""
    if len(a) == 0 or len(b) == 0:
        return np.nan

    u_stat = mannwhitneyu(
        a, b, alternative="two-sided", method="asymptotic"
    ).statistic

    return 2 * u_stat / (len(a) * len(b)) - 1


def bh_adjust(p_values):
    """Benjamini-Hochberg multiple-testing correction."""
    p = np.asarray(p_values, dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.where(np.isfinite(p))[0]

    if len(valid) == 0:
        return result

    values = p[valid]
    order = np.argsort(values)
    sorted_p = values[order]
    m = len(sorted_p)

    adjusted = sorted_p * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)

    result[valid[order]] = adjusted
    return result


def main():
    truth = pd.read_csv(GROUNDTRUTH)
    predictions = pd.read_csv(PREDICTIONS)

    keys = ["chr", "start", "end"]

    required_truth = set(keys + ["label", "young_mean"])
    required_predictions = set(keys + ["pred_atac_mean"])

    if not required_truth.issubset(truth.columns):
        raise ValueError(
            f"Missing ground-truth columns: "
            f"{sorted(required_truth - set(truth.columns))}"
        )

    if not required_predictions.issubset(predictions.columns):
        raise ValueError(
            f"Missing prediction columns: "
            f"{sorted(required_predictions - set(predictions.columns))}"
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

    for col in ["young_mean", "pred_atac_mean"]:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    data = data.replace([np.inf, -np.inf], np.nan)
    data = data.dropna(
        subset=["label", "young_mean", "pred_atac_mean"]
    )

    # 1. Ranking direction for each one-vs-rest task.
    # Negative is expected to have lower, rather than higher, scores.
    ranking_rows = []

    for label in CLASSES:
        y_true = (data["label"] == label).astype(int)
        scores = data["pred_atac_mean"]

        auc_forward = roc_auc_score(y_true, scores)
        auc_reverse = roc_auc_score(y_true, -scores)

        ranking_rows.append({
            "target_class": label,
            "roc_auc_forward": auc_forward,
            "roc_auc_reverse": auc_reverse,
            "better_direction": (
                "higher_score" if auc_forward >= auc_reverse
                else "lower_score"
            ),
        })

    ranking = pd.DataFrame(ranking_rows)
    ranking.to_csv(OUT / "class_ranking_direction.csv", index=False)

    # 2. Pairwise effect sizes for predicted ATAC signal.
    effect_rows = []

    for class_a, class_b in combinations(CLASSES, 2):
        a = data.loc[
            data["label"] == class_a, "pred_atac_mean"
        ].to_numpy()
        b = data.loc[
            data["label"] == class_b, "pred_atac_mean"
        ].to_numpy()

        test = mannwhitneyu(a, b, alternative="two-sided")

        effect_rows.append({
            "class_a": class_a,
            "class_b": class_b,
            "n_a": len(a),
            "n_b": len(b),
            "median_a": np.median(a),
            "median_b": np.median(b),
            "cliffs_delta_a_vs_b": cliff_delta(a, b),
            "p_value": test.pvalue,
        })

    effects = pd.DataFrame(effect_rows)
    effects["p_adjusted_BH"] = bh_adjust(effects["p_value"])
    effects.to_csv(OUT / "class_effect_sizes.csv", index=False)

    # 3. Exploratory baseline adjustment for RS_gained vs RS_lost.
    # Fit a pooled linear relationship: prediction ~ young_mean.
    # Compare residual predictions between gained and lost regions.
    # This is a sensitivity analysis, not a causal or state-aware model.
    x = data["young_mean"].to_numpy(dtype=float)
    y = data["pred_atac_mean"].to_numpy(dtype=float)

    design = np.column_stack([np.ones(len(data)), x])
    coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coefficients
    data["baseline_adjusted_prediction"] = y - fitted

    gained = data.loc[
        data["label"] == "RS_gained", "baseline_adjusted_prediction"
    ].to_numpy()
    lost = data.loc[
        data["label"] == "RS_lost", "baseline_adjusted_prediction"
    ].to_numpy()

    test = mannwhitneyu(gained, lost, alternative="two-sided")

    adjusted_result = pd.DataFrame([{
        "comparison": "RS_gained_vs_RS_lost",
        "adjustment": "pooled_linear_prediction_from_young_mean",
        "n_gained": len(gained),
        "n_lost": len(lost),
        "median_residual_gained": np.median(gained),
        "median_residual_lost": np.median(lost),
        "cliffs_delta_gained_vs_lost": cliff_delta(gained, lost),
        "mann_whitney_p_value": test.pvalue,
        "interpretation_note": (
            "Exploratory residual comparison; not a causal test and "
            "does not make the sequence model senescence-aware."
        ),
    }])

    adjusted_result.to_csv(
        OUT / "gained_lost_baseline_adjustment.csv", index=False
    )

    print("Follow-up evaluation completed.")
    print(f"Regions evaluated: {len(data)}")

    print("\nRanking direction:")
    print(ranking.to_string(index=False))

    print("\nPairwise effect sizes:")
    print(
        effects[
            [
                "class_a", "class_b", "cliffs_delta_a_vs_b",
                "p_adjusted_BH"
            ]
        ].to_string(index=False)
    )

    print("\nBaseline-adjusted RS_gained vs RS_lost:")
    print(adjusted_result.to_string(index=False))

    print(f"\nTables saved in: {OUT}")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()
