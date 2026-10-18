from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score, mean_absolute_error

ROOT = Path(__file__).resolve().parents[2]
TRUTH = ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
PRED = ROOT / "1_data_raw/3_alphagenome/02_atac_run/benchmark_10000_alphagenome_atac_pred.csv"
OUT = ROOT / "3_results/tables/change_association"
OUT.mkdir(parents=True, exist_ok=True)

keys = ["chr", "start", "end"]
truth = pd.read_csv(TRUTH)
pred = pd.read_csv(PRED)
truth["chr"] = truth["chr"].astype(str)
pred["chr"] = pred["chr"].astype(str)

df = truth.merge(
    pred[keys + ["pred_atac_mean"]],
    on=keys,
    how="inner",
    validate="one_to_one",
)

numeric = ["young_mean", "rs_mean", "log2FC_rs_vs_young", "pred_atac_mean"]
for col in numeric:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=numeric + ["label"])

# Exclude synthetic negatives: they have no experimentally measured
# young-to-senescent change suitable for this regression question.
df = df[df["label"].isin(["RS_gained", "RS_lost", "Stable"])].copy()

if df["chr"].nunique() < 5:
    raise ValueError("Need at least five chromosomes for five-fold grouped CV.")

y = df["log2FC_rs_vs_young"].to_numpy()
groups = df["chr"].to_numpy()

features = {
    "young_only": ["young_mean"],
    "young_plus_prediction": ["young_mean", "pred_atac_mean"],
}

splitter = GroupKFold(n_splits=5)
oof = {name: np.full(len(df), np.nan) for name in features}

for train_idx, test_idx in splitter.split(df, y, groups):
    for name, cols in features.items():
        model = LinearRegression()
        model.fit(df.iloc[train_idx][cols], y[train_idx])
        oof[name][test_idx] = model.predict(df.iloc[test_idx][cols])

rows = []
for name, predictions in oof.items():
    rows.append({
        "model": name,
        "n": len(df),
        "chromosomes": df["chr"].nunique(),
        "out_of_fold_R2": r2_score(y, predictions),
        "out_of_fold_MAE": mean_absolute_error(y, predictions),
    })

results = pd.DataFrame(rows)
results.to_csv(OUT / "chromosome_grouped_cv.csv", index=False)

print("Chromosome-grouped 5-fold cross-validation")
print("Regions:", len(df))
print("Chromosomes:", df["chr"].nunique())
print(results.round(5).to_string(index=False))

delta_r2 = (
    results.loc[results["model"] == "young_plus_prediction", "out_of_fold_R2"].iloc[0]
    - results.loc[results["model"] == "young_only", "out_of_fold_R2"].iloc[0]
)
delta_mae = (
    results.loc[results["model"] == "young_only", "out_of_fold_MAE"].iloc[0]
    - results.loc[results["model"] == "young_plus_prediction", "out_of_fold_MAE"].iloc[0]
)
print(f"\nChange in out-of-fold R2 (extended - baseline): {delta_r2:.5f}")
print(f"MAE improvement (baseline - extended): {delta_mae:.5f}")
print("Saved:", OUT / "chromosome_grouped_cv.csv")
print(
    "\nCaution: held-out chromosomes test transfer within this selected benchmark; "
    "they do not remove selection bias or replace validation on an independent dataset."
)
