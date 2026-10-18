from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm

root = Path(__file__).resolve().parents[2]
folder = root / "1_data_raw/5_motif_analysis"

df = pd.read_csv(folder / "tead1_motif_diagnostic_regions.csv")
df = df.loc[df["valid_fraction"] > 0].copy()

df["has_motif_match"] = df["has_motif_match"].astype(int)
df["log_length"] = np.log(df["region_length"].clip(lower=1))

for col in ["log_length", "gc_fraction"]:
    sd = df[col].std()
    if sd == 0:
        raise ValueError(f"No variation in {col}")
    df[col + "_z"] = (df[col] - df[col].mean()) / sd

# Stable is the reference class.
X = np.column_stack([
    np.ones(len(df)),
    (df["label"] == "Negative").astype(float),
    (df["label"] == "RS_gained").astype(float),
    (df["label"] == "RS_lost").astype(float),
    df["log_length_z"],
    df["gc_fraction_z"]
])
y = df["has_motif_match"].to_numpy(dtype=float)

def objective(beta):
    eta = np.clip(X @ beta, -30, 30)
    return np.sum(np.logaddexp(0, eta) - y * eta)

def gradient(beta):
    eta = np.clip(X @ beta, -30, 30)
    p = 1 / (1 + np.exp(-eta))
    return X.T @ (p - y)

fit = minimize(
    objective, np.zeros(X.shape[1]),
    jac=gradient, method="BFGS"
)

beta = fit.x
p = 1 / (1 + np.exp(-np.clip(X @ beta, -30, 30)))
information = X.T @ ((p * (1 - p))[:, None] * X)
cov = np.linalg.pinv(information)
se = np.sqrt(np.maximum(np.diag(cov), 0))

names = [
    "Intercept", "Negative_vs_Stable", "RS_gained_vs_Stable",
    "RS_lost_vs_Stable", "Log_length_per_SD", "GC_per_SD"
]
rows = []

for i, name in enumerate(names):
    z = beta[i] / se[i] if se[i] > 0 else np.nan
    rows.append({
        "term": name,
        "coefficient": beta[i],
        "standard_error": se[i],
        "odds_ratio": np.exp(beta[i]),
        "ci95_lower": np.exp(beta[i] - 1.96 * se[i]),
        "ci95_upper": np.exp(beta[i] + 1.96 * se[i]),
        "p_value": 2 * norm.sf(abs(z)) if np.isfinite(z) else np.nan
    })

results = pd.DataFrame(rows)
out = folder / "tead1_adjusted_logistic_results.csv"
results.to_csv(out, index=False)

print("Optimizer converged:", fit.success)
print("Optimizer message:", fit.message)
print("Evaluable regions:", len(df))
print("\nObserved motif matches by class:")
print(df.groupby("label")["has_motif_match"].agg(
    ["count", "sum", "mean"]
).round(4))
print("\nAdjusted logistic regression (Stable reference):")
print(results.round(5).to_string(index=False))
print("\nSaved:", out)
print("\nCaution: the motif threshold remains exploratory.")
