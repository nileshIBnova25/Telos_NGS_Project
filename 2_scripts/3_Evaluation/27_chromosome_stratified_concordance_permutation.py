import os
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

SEED = 20261024
N_PERMUTATIONS = 20000

base = "3_results/tables/integration/"
gene_file = base + "atac_rnaseq_single_class_gene_concordance_v2.csv"
association_file = base + "atac_rnaseq_direction_concordance_associations.csv"

genes = pd.read_csv(gene_file)
associations = pd.read_csv(association_file)

# Map each gene to its annotation chromosome.
chromosomes = associations[
    ["ensembl_id", "annotation_chr"]
].dropna().drop_duplicates()

if chromosomes["ensembl_id"].duplicated().any():
    raise ValueError("At least one gene maps to multiple chromosomes.")

genes = genes.merge(
    chromosomes, on="ensembl_id", how="left", validate="one_to_one"
)

if genes["annotation_chr"].isna().any():
    raise ValueError("Some genes have no annotation chromosome.")

genes = genes[
    genes["atac_class"].notna()
    & genes["concordant"].notna()
].copy()

genes["concordant"] = genes["concordant"].astype(bool)

classes = sorted(genes["atac_class"].unique())
if "RS_gained" not in classes:
    raise ValueError("RS_gained class is missing.")

rng = np.random.default_rng(SEED)

# Keep the observed class counts within each chromosome fixed.
groups = []
for chrom, group in genes.groupby("annotation_chr", sort=True):
    groups.append({
        "chromosome": chrom,
        "indices": group.index.to_numpy(),
        "labels": group["atac_class"].to_numpy().copy(),
        "outcomes": group["concordant"].to_numpy(dtype=bool)
    })

def rates_for_labels(labels, outcomes):
    result = {}
    for cls in classes:
        mask = labels == cls
        result[cls] = outcomes[mask].mean() if mask.sum() else np.nan
    return result

observed_rates = rates_for_labels(
    genes["atac_class"].to_numpy(),
    genes["concordant"].to_numpy(dtype=bool)
)

def gained_minus_class(rates, cls):
    return rates["RS_gained"] - rates[cls]

comparisons = [c for c in classes if c != "RS_gained"]
observed_diffs = {
    cls: gained_minus_class(observed_rates, cls)
    for cls in comparisons
}

# Also compare gained genes with all non-gained genes pooled.
observed_labels = genes["atac_class"].to_numpy()
observed_outcomes = genes["concordant"].to_numpy(dtype=bool)
gained_mask = observed_labels == "RS_gained"

observed_pooled_diff = (
    observed_outcomes[gained_mask].mean()
    - observed_outcomes[~gained_mask].mean()
)

exceed_pair = {cls: 0 for cls in comparisons}
exceed_pooled = 0

for _ in range(N_PERMUTATIONS):
    permuted_labels = observed_labels.copy()

    for group in groups:
        permuted_labels[group["indices"]] = rng.permutation(group["labels"])

    rates = rates_for_labels(permuted_labels, observed_outcomes)

    for cls in comparisons:
        diff = gained_minus_class(rates, cls)
        if diff >= observed_diffs[cls]:
            exceed_pair[cls] += 1

    perm_gained = permuted_labels == "RS_gained"
    pooled_diff = (
        observed_outcomes[perm_gained].mean()
        - observed_outcomes[~perm_gained].mean()
    )
    if pooled_diff >= observed_pooled_diff:
        exceed_pooled += 1

# Add-one correction avoids reporting an empirical p-value of exactly zero.
rows = []
for cls in comparisons:
    rows.append({
        "analysis": f"RS_gained_vs_{cls}",
        "gained_rate": observed_rates["RS_gained"],
        "comparison_rate": observed_rates[cls],
        "observed_difference": observed_diffs[cls],
        "permutations": N_PERMUTATIONS,
        "empirical_one_sided_p": (
            exceed_pair[cls] + 1
        ) / (N_PERMUTATIONS + 1)
    })

rows.append({
    "analysis": "RS_gained_vs_all_other_classes",
    "gained_rate": observed_outcomes[gained_mask].mean(),
    "comparison_rate": observed_outcomes[~gained_mask].mean(),
    "observed_difference": observed_pooled_diff,
    "permutations": N_PERMUTATIONS,
    "empirical_one_sided_p": (
        exceed_pooled + 1
    ) / (N_PERMUTATIONS + 1)
})

results = pd.DataFrame(rows)
results["p_adjusted_BH"] = multipletests(
    results["empirical_one_sided_p"], method="fdr_bh"
)[1]

out_file = base + "atac_rnaseq_chromosome_stratified_permutation.csv"
results.to_csv(out_file, index=False)

print("===== ANALYSIS DESIGN =====")
print("Directionally classifiable single-class genes:", len(genes))
print("Chromosomes:", genes["annotation_chr"].nunique())
print("Permutations:", N_PERMUTATIONS)
print("Random seed:", SEED)
print("Class counts preserved within each chromosome.")
print("Gene concordance outcomes remain fixed.")
print("\n===== OBSERVED CONCORDANCE RATES =====")
for cls in classes:
    n = int((genes["atac_class"] == cls).sum())
    print(f"{cls}: {observed_rates[cls] * 100:.2f}% (n={n})")

print("\n===== EMPIRICAL PERMUTATION TESTS =====")
print(results.to_string(index=False))
print("\nSaved:", out_file)
print(
    "\nCaution: chromosome stratification does not preserve local "
    "spatial dependence or all biological relationships among genes."
)
