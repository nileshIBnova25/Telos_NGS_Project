import os
import itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, permutation_test
from statsmodels.stats.multitest import multipletests

INPUT = "3_results/tables/unified_motif_analysis/all_nine_motifs_combined.csv"
OUTDIR = "3_results/tables/motif_cooccurrence_score_correlations/class_tests"
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(INPUT)

# Use the calibrated motif score excess for each region and motif.
wide = df.pivot(index="region_id", columns="motif", values="score_excess")
labels = df.drop_duplicates("region_id").set_index("region_id")["label"]

wide["label"] = labels
motifs = [c for c in wide.columns if c != "label"]
classes = ["RS_gained", "RS_lost", "Stable", "Negative"]

print("Regions with complete motif scores:", wide[motifs].dropna().shape[0])
print("Motifs:", len(motifs))
print("Classes:", wide["label"].value_counts().to_dict())

def correlation_difference(a, b, c, d):
    # Fisher transformation for comparing independent correlations.
    r1 = np.clip(spearmanr(a, b).statistic, -0.999999, 0.999999)
    r2 = np.clip(spearmanr(c, d).statistic, -0.999999, 0.999999)
    return np.arctanh(r1) - np.arctanh(r2), r1, r2

rows = []
for m1, m2 in itertools.combinations(motifs, 2):
    for c1, c2 in itertools.combinations(classes, 2):
        a = wide.loc[wide.label == c1, m1].dropna()
        # Keep paired values for each class.
        d1 = wide.loc[wide.label == c1, [m1, m2]].dropna()
        d2 = wide.loc[wide.label == c2, [m1, m2]].dropna()

        if len(d1) < 4 or len(d2) < 4:
            continue

        z_diff, r1, r2 = correlation_difference(
            d1[m1], d1[m2], d2[m1], d2[m2]
        )
        se = np.sqrt(1/(len(d1)-3) + 1/(len(d2)-3))
        z = z_diff / se
        # Two-sided normal approximation.
        from scipy.stats import norm
        p = 2 * norm.sf(abs(z))

        rows.append({
            "motif_1": m1,
            "motif_2": m2,
            "class_1": c1,
            "class_2": c2,
            "n_class_1": len(d1),
            "n_class_2": len(d2),
            "rho_class_1": r1,
            "rho_class_2": r2,
            "rho_difference": r1-r2,
            "p_value": p
        })

res = pd.DataFrame(rows)
if len(res):
    res["BH_q_value"] = multipletests(res["p_value"], method="fdr_bh")[1]
    res = res.sort_values("BH_q_value")
    out = os.path.join(OUTDIR, "motif_pair_correlation_class_differences.csv")
    res.to_csv(out, index=False)
    print("\nSaved:", out)
    print("Tests:", len(res))
    print("FDR-significant tests (q < 0.05):", (res.BH_q_value < 0.05).sum())
    print("\nTop 15 results:")
    print(res.head(15).to_string(index=False))
else:
    print("No valid class comparisons were produced.")
