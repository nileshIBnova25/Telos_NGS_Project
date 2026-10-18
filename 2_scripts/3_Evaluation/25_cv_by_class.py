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

cols = ["young_mean", "log2FC_rs_vs_young", "pred_atac_mean"]
for col in cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=cols + ["label"])
df = df[df["label"].isin(["RS_gained", "RS_lost", "Stable"])].copy()
df = df.reset_index(drop=True)

y = df["log2FC_rs_vs_young"].to_numpy()
groups = df["chr"].to_numpy()
features = {
    "young_only": ["young_mean"],
    "young_plus_prediction": ["young_mean", "pred_atac_mean"],
}
oof = {name: np.full(len(df), np.nan) for name in features}

cv = GroupKFold(n_splits=5)
for train_idx, test_idx in cv.split(df, y, groups):
    for name, columns in features.items():
        model = LinearRegression()
        model.fit(df.iloc[train_idx][columns], y[train_idx])
        oof[name][test_idx] = model.predict(df.iloc[test_idx][columns])

df["oof_young_only"] = oof["young_only"]
df["oof_young_plus_prediction"] = oof["young_plus_prediction"]
df.to_csv(OUT / "chromosome_grouped_cv_predictions.csv", index=False)

rows = []
for label, group in df.groupby("label"):
    actual = group["log2FC_rs_vs_young"].to_numpy()
    for model_name, pred_col in [
        ("young_only", "oof_young_only"),
        ("young_plus_prediction", "oof_young_plus_prediction"),
    ]:
        predicted = group[pred_col].to_numpy()
        rows.append({
            "class": label,
            "model": model_name,
            "n": len(group),
            "out_of_fold_R2": r2_score(actual, predicted),
            "out_of_fold_MAE": mean_absolute_error(actual, predicted),
        })

results = pd.DataFrame(rows)
results.to_csv(OUT / "chromosome_grouped_cv_by_class.csv", index=False)

print("Chromosome-grouped cross-validation by class")
print(results.round(5).to_string(index=False))

print("\nExtended-model improvement by class (positive means improvement):")
for label in ["RS_gained", "RS_lost", "Stable"]:
    subset = results[results["class"] == label].set_index("model")
    delta_r2 = (
        subset.loc["young_plus_prediction", "out_of_fold_R2"]
        - subset.loc["young_only", "out_of_fold_R2"]
    )
    delta_mae = (
        subset.loc["young_only", "out_of_fold_MAE"]
        - subset.loc["young_plus_prediction", "out_of_fold_MAE"]
    )
    print(f"{label}: delta R2={delta_r2:.5f}; MAE improvement={delta_mae:.5f}")

print("\nSaved per-region predictions and class metrics to:", OUT)
print(
    "Caution: classes were selected using experimental fold changes; "
    "class-wise results are exploratory and not independent genome-wide validation."
)
