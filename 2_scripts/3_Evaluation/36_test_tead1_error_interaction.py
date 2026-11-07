from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

root = Path.cwd()

data_dir = root / "1_data_raw/5_motif_analysis"
pred_path = root / "1_data_raw/3_alphagenome/02_atac_run/benchmark_10000_alphagenome_atac_pred.csv"
truth_path = root / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
cal_path = data_dir / "tead1_dinucleotide_calibration_regions.csv"
diag_path = data_dir / "tead1_motif_diagnostic_regions.csv"

out_dir = root / "3_results/tables/tead1_alphagenome_error"
out_dir.mkdir(parents=True, exist_ok=True)

pred = pd.read_csv(pred_path)
truth = pd.read_csv(truth_path)
cal = pd.read_csv(cal_path)
diag = pd.read_csv(diag_path)

cal["tead1_score_excess"] = (
    cal["observed_max_score"] - cal["null_mean_max_score"]
)

df = pred.merge(
    truth[["name", "rs_mean"]],
    on="name", how="inner", validate="one_to_one"
)
df = df.merge(
    cal[["region_id", "tead1_score_excess"]],
    left_on="name", right_on="region_id",
    how="inner", validate="one_to_one"
)
df = df.merge(
    diag[["region_id", "region_length", "gc_fraction", "valid_fraction"]],
    on="region_id", how="inner", validate="one_to_one"
)

df = df[
    (df["valid_fraction"] > 0) &
    df["tead1_score_excess"].notna() &
    df["gc_fraction"].notna() &
    df["region_length"].notna()
].copy()

df["signed_residual"] = df["pred_atac_mean"] - df["rs_mean"]
df["absolute_error"] = df["signed_residual"].abs()
df["log_length"] = np.log(df["region_length"].clip(lower=1))

# Standardize continuous predictors for comparable model coefficients.
for col in ["tead1_score_excess", "log_length", "gc_fraction"]:
    sd = df[col].std()
    if not np.isfinite(sd) or sd == 0:
        raise ValueError(f"Cannot standardize {col}")
    df[col + "_z"] = (df[col] - df[col].mean()) / sd

df["label"] = pd.Categorical(
    df["label"],
    categories=["Stable", "Negative", "RS_gained", "RS_lost"]
)

print("Matched, evaluable regions:", len(df))
print("Class counts:")
print(df["label"].value_counts().to_string())

# Save the exact analysis cohort for reproducibility.
df.to_csv(out_dir / "tead1_interaction_analysis_cohort.csv", index=False)

joint_rows = []
slope_rows = []

for outcome in ["absolute_error", "signed_residual"]:
    formula = (
        f"{outcome} ~ C(label) * tead1_score_excess_z"
        " + log_length_z + gc_fraction_z"
    )
    model = smf.ols(formula, data=df).fit(cov_type="HC3")

    # Joint robust Wald test: are all class-by-score interactions zero?
    interaction_names = [
        name for name in model.params.index
        if ":" in name and "tead1_score_excess_z" in name
    ]

    R = np.zeros((len(interaction_names), len(model.params)))
    for i, name in enumerate(interaction_names):
        R[i, list(model.params.index).index(name)] = 1

    test = model.wald_test(R, scalar=True)
    joint_rows.append({
        "outcome": outcome,
        "n": int(model.nobs),
        "interaction_terms": len(interaction_names),
        "robust_wald_chi2": float(test.statistic),
        "df": len(interaction_names),
        "p_value": float(test.pvalue)
    })

    print(f"\n========== {outcome} ==========")
    print("Formula:", formula)
    print("R-squared:", round(model.rsquared, 5))
    print("\nRobust regression coefficients:")
    print(model.summary().tables[1])

    # Estimate the TEAD1-score slope within each class.
    for label in ["Stable", "Negative", "RS_gained", "RS_lost"]:
        terms = ["tead1_score_excess_z"]
        if label != "Stable":
            terms.append(
                f"C(label)[T.{label}]:tead1_score_excess_z"
            )

        contrast = np.zeros(len(model.params))
        for term in terms:
            if term not in model.params.index:
                raise ValueError(f"Missing expected coefficient: {term}")
            contrast[list(model.params.index).index(term)] = 1

        t = model.t_test(contrast)
        slope_rows.append({
            "outcome": outcome,
            "label": label,
            "n": int((df["label"] == label).sum()),
            "slope_per_SD_score": float(np.asarray(t.effect).reshape(-1)[0]),
            "SE_HC3": float(np.asarray(t.sd).reshape(-1)[0]),
            "p_value": float(np.asarray(t.pvalue).reshape(-1)[0])
        })

joint = pd.DataFrame(joint_rows)
joint["BH_q_value"] = multipletests(joint["p_value"], method="fdr_bh")[1]
joint.to_csv(out_dir / "tead1_error_interaction_joint_tests.csv", index=False)

slopes = pd.DataFrame(slope_rows)
slopes["BH_q_value"] = multipletests(slopes["p_value"], method="fdr_bh")[1]
slopes.to_csv(out_dir / "tead1_error_interaction_class_slopes.csv", index=False)

print("\n========== JOINT CLASS-BY-SCORE INTERACTION TESTS ==========")
print(joint.round(6).to_string(index=False))

print("\n========== SCORE SLOPE WITHIN EACH CLASS ==========")
print(slopes.round(6).to_string(index=False))

print("\nSaved results to:", out_dir)
print(
    "\nInterpretation note: HC3 robust standard errors are used. "
    "The joint test evaluates whether the score-error slopes differ "
    "across classes, controlling for log region length and GC fraction. "
    "Associations do not establish YAP1 binding or causality."
)
