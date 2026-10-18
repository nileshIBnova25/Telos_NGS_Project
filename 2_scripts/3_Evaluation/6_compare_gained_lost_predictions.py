from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]

TRUTH_FILE = (
    ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)
PRED_FILE = (
    ROOT
    / "1_data_raw/3_alphagenome/02_atac_run/"
      "benchmark_10000_alphagenome_atac_pred.csv"
)
OUT_FILE = ROOT / "3_results/tables/gained_vs_lost_prediction_comparison.csv"

RANDOM_SEED = 42
N_BOOTSTRAP = 2000

truth = pd.read_csv(TRUTH_FILE)
pred = pd.read_csv(PRED_FILE)

df = truth.merge(
    pred[["chr", "start", "end", "pred_atac_mean"]],
    on=["chr", "start", "end"],
    how="inner",
    validate="one_to_one",
)

df = df[df["label"].isin(["RS_gained", "RS_lost"])].copy()
df = df.dropna(subset=["pred_atac_mean"])

gained = df.loc[
    df["label"] == "RS_gained", "pred_atac_mean"
].to_numpy()

lost = df.loc[
    df["label"] == "RS_lost", "pred_atac_mean"
].to_numpy()

if len(gained) == 0 or len(lost) == 0:
    raise ValueError("Both RS_gained and RS_lost must contain observations.")

# Positive class: RS_gained. A higher score is treated as evidence
# for RS_gained, and a lower score as evidence for RS_lost.
scores = np.concatenate([gained, lost])
labels = np.concatenate([
    np.ones(len(gained), dtype=int),
    np.zeros(len(lost), dtype=int),
])

auc = roc_auc_score(labels, scores)

# Bootstrap the two classes independently to estimate uncertainty.
rng = np.random.default_rng(RANDOM_SEED)
bootstrap_aucs = []

for _ in range(N_BOOTSTRAP):
    sampled_gained = rng.choice(gained, size=len(gained), replace=True)
    sampled_lost = rng.choice(lost, size=len(lost), replace=True)

    sampled_scores = np.concatenate([sampled_gained, sampled_lost])
    sampled_labels = np.concatenate([
        np.ones(len(sampled_gained), dtype=int),
        np.zeros(len(sampled_lost), dtype=int),
    ])

    bootstrap_aucs.append(
        roc_auc_score(sampled_labels, sampled_scores)
    )

ci_low, ci_high = np.percentile(bootstrap_aucs, [2.5, 97.5])

# Mann-Whitney test.
u_stat, p_value = mannwhitneyu(
    gained, lost, alternative="two-sided"
)

# Cliff's delta: P(gained > lost) - P(gained < lost).
# The AUC also includes half credit for ties.
cliffs_delta = 2 * auc - 1

results = pd.DataFrame([{
    "class_1": "RS_gained",
    "class_2": "RS_lost",
    "n_gained": len(gained),
    "n_lost": len(lost),
    "median_pred_gained": np.median(gained),
    "median_pred_lost": np.median(lost),
    "mean_pred_gained": np.mean(gained),
    "mean_pred_lost": np.mean(lost),
    "roc_auc_gained_vs_lost": auc,
    "bootstrap_95ci_lower": ci_low,
    "bootstrap_95ci_upper": ci_high,
    "mann_whitney_u": u_stat,
    "mann_whitney_p": p_value,
    "cliffs_delta": cliffs_delta,
    "bootstrap_replicates": N_BOOTSTRAP,
    "random_seed": RANDOM_SEED,
}])

OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(OUT_FILE, index=False)

print("RS_GAINED VS RS_LOST PREDICTION COMPARISON")
print("-----------------------------------------")
print(f"RS_gained regions: {len(gained)}")
print(f"RS_lost regions:   {len(lost)}")
print(f"Median prediction, RS_gained: {np.median(gained):.6f}")
print(f"Median prediction, RS_lost:   {np.median(lost):.6f}")
print(f"ROC-AUC: {auc:.6f}")
print(f"Bootstrap 95% CI: [{ci_low:.6f}, {ci_high:.6f}]")
print(f"Mann-Whitney U: {u_stat:.1f}")
print(f"Mann-Whitney p-value: {p_value:.6e}")
print(f"Cliff's delta: {cliffs_delta:.6f}")
print(f"\nSaved results to: {OUT_FILE}")
