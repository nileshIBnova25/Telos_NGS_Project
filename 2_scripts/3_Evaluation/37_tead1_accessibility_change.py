from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

root = Path.cwd()

motif_path = root / "1_data_raw/5_motif_analysis/tead1_motif_region_scores.csv"
truth_path = root / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"

out_dir = root / "3_results/tables/tead1_accessibility_change"
out_dir.mkdir(parents=True, exist_ok=True)

motif = pd.read_csv(motif_path)
truth = pd.read_csv(truth_path)

df = motif.merge(
    truth[["name", "log2FC_rs_vs_young", "young_mean", "rs_mean"]],
    left_on="region_id",
    right_on="name",
    how="inner",
    validate="one_to_one"
)

df = df.dropna(subset=["log2FC_rs_vs_young", "has_motif_match"]).copy()
df["has_motif_match"] = df["has_motif_match"].astype(int)

print("Matched regions:", len(df))
print("\nRegion counts by class and motif status:")
print(pd.crosstab(df["label"], df["has_motif_match"]).to_string())

# Compare observed accessibility changes for motif-present vs motif-absent regions.
rows = []

for label, group in df.groupby("label"):
    present = group.loc[group["has_motif_match"] == 1, "log2FC_rs_vs_young"]
    absent = group.loc[group["has_motif_match"] == 0, "log2FC_rs_vs_young"]

    if len(present) == 0 or len(absent) == 0:
        continue

    test = mannwhitneyu(present, absent, alternative="two-sided")

    rows.append({
        "label": label,
        "n_motif_present": len(present),
        "n_motif_absent": len(absent),
        "median_log2FC_present": present.median(),
        "median_log2FC_absent": absent.median(),
        "mean_log2FC_present": present.mean(),
        "mean_log2FC_absent": absent.mean(),
        "median_difference_present_minus_absent": (
            present.median() - absent.median()
        ),
        "mann_whitney_p": test.pvalue
    })

results = pd.DataFrame(rows)
results["BH_q_value"] = multipletests(
    results["mann_whitney_p"], method="fdr_bh"
)[1]

results.to_csv(out_dir / "tead1_match_accessibility_by_class.csv", index=False)

print("\nAccessibility change by motif status within each class:")
print(results.round(6).to_string(index=False))

# Overall comparison
present = df.loc[df["has_motif_match"] == 1, "log2FC_rs_vs_young"]
absent = df.loc[df["has_motif_match"] == 0, "log2FC_rs_vs_young"]
overall = mannwhitneyu(present, absent, alternative="two-sided")

print("\nOverall comparison:")
print("Motif present:", len(present))
print("Motif absent:", len(absent))
print("Median log2FC, present:", round(present.median(), 6))
print("Median log2FC, absent:", round(absent.median(), 6))
print("Median difference:", round(present.median() - absent.median(), 6))
print("Mann-Whitney p-value:", overall.pvalue)

# Continuous score association, overall and by class
corr_rows = []

for label, group in [("ALL", df)] + list(df.groupby("label")):
    rho, p = spearmanr(
        group["max_motif_score"],
        group["log2FC_rs_vs_young"]
    )
    corr_rows.append({
        "label": label,
        "n": len(group),
        "spearman_rho": rho,
        "p_value": p
    })

corr = pd.DataFrame(corr_rows)
corr["BH_q_value"] = multipletests(corr["p_value"], method="fdr_bh")[1]
corr.to_csv(out_dir / "tead1_score_accessibility_correlations.csv", index=False)

print("\nContinuous motif score vs observed log2FC:")
print(corr.round(6).to_string(index=False))

print("\nSaved results to:", out_dir)
