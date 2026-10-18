from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr, spearmanr

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
print("Matched usable regions:", len(df))

# 1. Correlation with the experimentally observed log2 fold change.
corr_rows = []
groups = [("ALL", df)] + [
    (label, group)
    for label, group in df.groupby("label")
    if label != "Negative"
]

for label, group in groups:
    if len(group) < 3:
        continue
    corr_rows.append({
        "group": label,
        "n": len(group),
        "pearson_r_pred_vs_log2FC": pearsonr(
            group["pred_atac_mean"], group["log2FC_rs_vs_young"]
        ).statistic,
        "spearman_r_pred_vs_log2FC": spearmanr(
            group["pred_atac_mean"], group["log2FC_rs_vs_young"]
        ).statistic,
    })

corr = pd.DataFrame(corr_rows)
corr.to_csv(OUT / "prediction_vs_change_correlations.csv", index=False)

# 2. Exploratory regression:
# Does predicted signal associate with log2 fold change after accounting
# for young signal? Standardize predictors for coefficient comparability.
# This is not causal and is sensitive to benchmark selection.
model_data = df[df["label"] != "Negative"].copy()

for col in ["pred_atac_mean", "young_mean"]:
    sd = model_data[col].std(ddof=0)
    if sd == 0:
        raise ValueError(f"No variation in {col}")
    model_data[col + "_z"] = (
        model_data[col] - model_data[col].mean()
    ) / sd

model_data["log2FC_rs_vs_young"] = pd.to_numeric(
    model_data["log2FC_rs_vs_young"], errors="coerce"
)

# Use HC3 robust standard errors. No claim of independent observations is made.
X_base = sm.add_constant(model_data[["young_mean_z"]])
X_full = sm.add_constant(model_data[["young_mean_z", "pred_atac_mean_z"]])
y = model_data["log2FC_rs_vs_young"]

base = sm.OLS(y, X_base).fit(cov_type="HC3")
full = sm.OLS(y, X_full).fit(cov_type="HC3")

reg_rows = []
for model_name, fitted in [("young_only", base), ("young_plus_prediction", full)]:
    for term in fitted.params.index:
        reg_rows.append({
            "model": model_name,
            "term": term,
            "coefficient": fitted.params[term],
            "std_error_HC3": fitted.bse[term],
            "p_value_HC3": fitted.pvalues[term],
            "r_squared": fitted.rsquared,
            "n": int(fitted.nobs),
        })

reg = pd.DataFrame(reg_rows)
reg.to_csv(OUT / "exploratory_regression.csv", index=False)

print("\nCorrelation with observed log2 fold change:")
print(corr.round(4).to_string(index=False))
print("\nExploratory regression (HC3 standard errors):")
print(reg.round(6).to_string(index=False))
print("\nSaved results to:", OUT)
print(
    "\nCaution: classes were selected using experimental fold changes; "
    "the overall association and regression are exploratory, not independent validation."
)
