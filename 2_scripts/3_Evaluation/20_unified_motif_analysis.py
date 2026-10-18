from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
CAL = ROOT / "3_results/tables/paper_guided_calibration"
TEAD = ROOT / "3_results/tables/tead1_calibration_1000"
DIAG = ROOT / "1_data_raw/5_motif_analysis/tead1_motif_diagnostic_regions.csv"
OUT = ROOT / "3_results/tables/unified_motif_analysis"
OUT.mkdir(parents=True, exist_ok=True)

# Read and standardize the nine motif calibration tables.
tead = pd.read_csv(
    TEAD / "tead1_dinucleotide_calibration_regions_1000.csv"
)
tead["motif"] = "TEAD1"

four = pd.read_csv(CAL / "four_motif_region_calibration.csv")
remaining = pd.read_csv(
    CAL / "remaining_four_motif_region_calibration.csv"
)

data = pd.concat([tead, four, remaining], ignore_index=True)

expected = {
    "TEAD1", "CEBPA", "FOS", "RELA", "SMAD4",
    "E2F1", "JUN", "NFKB1", "SMAD3"
}
if set(data["motif"].unique()) != expected:
    raise ValueError("The combined tables do not contain exactly the expected nine motifs.")

if data.duplicated(["motif", "region_id"]).any():
    raise ValueError("Duplicate motif-region records detected.")

# Join the existing sequence diagnostics without changing source files.
diag = pd.read_csv(DIAG)
diag = diag[["region_id", "region_length", "gc_fraction", "valid_fraction"]]

if diag["region_id"].duplicated().any():
    raise ValueError("Duplicate region IDs in diagnostics.")

data = data.merge(diag, on="region_id", how="left", validate="many_to_one")

if data["region_length"].isna().any():
    raise ValueError("Some motif records did not match sequence diagnostics.")

# Save the complete combined table, including unscorable records.
data["score_excess"] = (
    data["observed_max_score"] - data["null_mean_max_score"]
)
data.to_csv(OUT / "all_nine_motifs_combined.csv", index=False)

# Model RS_gained and RS_lost against Stable separately for each motif.
# Use scorable regions in these three classes only.
results = []

for motif in sorted(expected):
    subset = data[
        (data["motif"] == motif)
        & data["label"].isin(["RS_gained", "RS_lost", "Stable"])
    ].copy()

    subset = subset.dropna(
        subset=["score_excess", "region_length", "gc_fraction"]
    )

    subset = subset[
        (subset["region_length"] > 0)
        & subset["gc_fraction"].between(0, 1)
    ].copy()

    subset["log_length"] = np.log(subset["region_length"])

    # Standardize covariates within this motif's analysis subset.
    for col in ["log_length", "gc_fraction"]:
        sd = subset[col].std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"Cannot standardize {col} for {motif}.")
        subset[col + "_z"] = (subset[col] - subset[col].mean()) / sd

    subset["RS_gained"] = (subset["label"] == "RS_gained").astype(int)
    subset["RS_lost"] = (subset["label"] == "RS_lost").astype(int)

    X = subset[
        ["RS_gained", "RS_lost", "log_length_z", "gc_fraction_z"]
    ]
    X = sm.add_constant(X, has_constant="add")
    y = subset["score_excess"]

    model = sm.OLS(y, X).fit(cov_type="HC3")

    for contrast in ["RS_gained", "RS_lost"]:
        beta = float(model.params[contrast])
        se = float(model.bse[contrast])
        p = float(model.pvalues[contrast])
        ci_low = beta - 1.96 * se
        ci_high = beta + 1.96 * se

        results.append({
            "motif": motif,
            "contrast": f"{contrast}_vs_Stable",
            "n_regions": len(subset),
            "beta_score_excess": beta,
            "ci95_low": ci_low,
            "ci95_high": ci_high,
            "p_value": p,
            "model_r_squared": float(model.rsquared),
        })

results = pd.DataFrame(results)

# Correct all 18 planned tests together.
results["BH_q_value_18_tests"] = multipletests(
    results["p_value"].values, method="fdr_bh"
)[1]

results = results.sort_values(["motif", "contrast"])
results.to_csv(OUT / "adjusted_score_excess_tests.csv", index=False)

# Descriptive summary of observed-minus-null scores by motif and class.
summary = (
    data[data["label"].isin(["Negative", "RS_gained", "RS_lost", "Stable"])]
    .groupby(["motif", "label"], dropna=False)
    .agg(
        n_records=("region_id", "size"),
        n_scorable=("score_excess", "count"),
        mean_score_excess=("score_excess", "mean"),
        median_score_excess=("score_excess", "median"),
        mean_observed_score=("observed_max_score", "mean"),
        mean_null_score=("null_mean_max_score", "mean"),
    )
    .reset_index()
)
summary.to_csv(OUT / "score_excess_class_summary.csv", index=False)

print("\nAdjusted score-excess models:")
print(results.round(5).to_string(index=False))

print("\nClass summary:")
print(summary.round(4).to_string(index=False))

print("\nSaved outputs to:", OUT)
