from pathlib import Path
from itertools import combinations

import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "1_data_raw/5_motif_analysis"

INPUT = FOLDER / "tead1_dinucleotide_calibration_regions.csv"
OUTPUT = FOLDER / "tead1_class_difference_tests.csv"

df = pd.read_csv(INPUT)
df = df.dropna(
    subset=["observed_max_score", "null_mean_max_score"]
).copy()

df["score_excess"] = (
    df["observed_max_score"] - df["null_mean_max_score"]
)

classes = ["Negative", "RS_gained", "RS_lost", "Stable"]
groups = {
    label: df.loc[df["label"] == label, "score_excess"].to_numpy()
    for label in classes
}

h, p_global = kruskal(*(groups[label] for label in classes))

print("Kruskal-Wallis overall test")
print(f"H statistic: {h:.5f}")
print(f"P-value: {p_global:.6g}")

rows = []

for a, b in combinations(classes, 2):
    x, y = groups[a], groups[b]
    result = mannwhitneyu(x, y, alternative="two-sided")

    # Rank-biserial correlation: positive means values tend to be
    # higher in class a than class b.
    effect = 2 * result.statistic / (len(x) * len(y)) - 1

    rows.append({
        "class_a": a,
        "class_b": b,
        "n_a": len(x),
        "n_b": len(y),
        "median_a": float(pd.Series(x).median()),
        "median_b": float(pd.Series(y).median()),
        "rank_biserial_a_vs_b": effect,
        "p_value": result.pvalue,
    })

results = pd.DataFrame(rows)

# Benjamini-Hochberg adjustment across the six pairwise comparisons.
m = len(results)
order = results["p_value"].argsort().to_numpy()
sorted_p = results.loc[order, "p_value"].to_numpy()

adjusted = (sorted_p * m / range(1, m + 1))
adjusted = pd.Series(adjusted).iloc[::-1].cummin().iloc[::-1]
adjusted = adjusted.clip(upper=1).to_numpy()

q_values = [None] * m
for rank, index in enumerate(order):
    q_values[index] = adjusted[rank]

results["BH_q_value"] = q_values
results.to_csv(OUTPUT, index=False)

print("\nPairwise Mann-Whitney tests")
print(results.round(6).to_string(index=False))
print("\nSaved:", OUTPUT)
print(
    "\nInterpretation caution: these tests compare observed-minus-null "
    "scores. They do not prove TEAD1 binding, and the null scores come "
    "from only 20 shuffles per region."
)
