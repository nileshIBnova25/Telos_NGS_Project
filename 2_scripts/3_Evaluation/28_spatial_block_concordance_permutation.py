import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

SEED = 20261024
N_PERMUTATIONS = 20000

base = "3_results/tables/integration/"

gene_path = base + "atac_rnaseq_single_class_gene_concordance_v2.csv"
assoc_path = base + "atac_rnaseq_direction_concordance_associations.csv"

genes = pd.read_csv(gene_path)
assoc = pd.read_csv(assoc_path)

# Retrieve one genomic coordinate per gene.
coords = assoc[
    ["ensembl_id", "annotation_chr", "annotation_start", "annotation_end"]
].dropna().drop_duplicates()

# Aggregate multiple annotation intervals for a gene into its full genomic span.
# All intervals must be on the same chromosome.
chrom_counts = coords.groupby("ensembl_id")["annotation_chr"].nunique()
if (chrom_counts > 1).any():
    raise ValueError("Some genes map to multiple chromosomes.")

coords = coords.groupby("ensembl_id", as_index=False).agg(
    annotation_chr=("annotation_chr", "first"),
    annotation_start=("annotation_start", "min"),
    annotation_end=("annotation_end", "max")
)

df = genes.merge(
    coords,
    on="ensembl_id",
    how="left",
    validate="one_to_one"
)

if df[["annotation_chr", "annotation_start", "annotation_end"]].isna().any().any():
    raise ValueError("Some genes lack genomic coordinates.")

df = df.reset_index(drop=True)
df["concordant"] = df["concordant"].astype(bool)
df["midpoint"] = (df["annotation_start"] + df["annotation_end"]) / 2

classes = sorted(df["atac_class"].unique())
if "RS_gained" not in classes:
    raise ValueError("RS_gained class missing.")

comparisons = [c for c in classes if c != "RS_gained"]
observed_labels = df["atac_class"].to_numpy().copy()
outcomes = df["concordant"].to_numpy(dtype=bool)
gained = observed_labels == "RS_gained"

observed_rates = {
    c: outcomes[observed_labels == c].mean()
    for c in classes
}

observed_differences = {
    c: observed_rates["RS_gained"] - observed_rates[c]
    for c in comparisons
}

observed_pooled_difference = (
    outcomes[gained].mean() - outcomes[~gained].mean()
)

print("Genes:", len(df))
print("Observed class concordance:")
for c in classes:
    print(f"  {c}: {observed_rates[c] * 100:.2f}% "
          f"(n={(observed_labels == c).sum()})")

all_results = []

for block_size in [5_000_000, 10_000_000]:
    work = df.copy()
    work["block_number"] = (
        work["midpoint"] // block_size
    ).astype(int)
    work["block_id"] = (
        work["annotation_chr"].astype(str)
        + "_"
        + work["block_number"].astype(str)
    )

    block_indices = [
        group.index.to_numpy()
        for _, group in work.groupby("block_id", sort=True)
    ]

    sizes = work.groupby("block_id").size()

    rng = np.random.default_rng(SEED + block_size)
    exceed = {c: 0 for c in comparisons}
    exceed_pooled = 0

    for _ in range(N_PERMUTATIONS):
        perm_labels = observed_labels.copy()

        # Shuffle class labels only within each chromosome-specific block.
        for idx in block_indices:
            perm_labels[idx] = rng.permutation(observed_labels[idx])

        perm_gained = perm_labels == "RS_gained"
        gained_rate = outcomes[perm_gained].mean()

        for c in comparisons:
            mask = perm_labels == c
            if mask.any():
                difference = gained_rate - outcomes[mask].mean()
                if difference >= observed_differences[c]:
                    exceed[c] += 1

        pooled_difference = (
            gained_rate - outcomes[~perm_gained].mean()
        )
        if pooled_difference >= observed_pooled_difference:
            exceed_pooled += 1

    for c in comparisons:
        all_results.append({
            "block_size_mb": block_size // 1_000_000,
            "comparison": f"RS_gained_vs_{c}",
            "gained_concordance": observed_rates["RS_gained"],
            "comparison_concordance": observed_rates[c],
            "observed_difference": observed_differences[c],
            "permutations": N_PERMUTATIONS,
            "empirical_one_sided_p": (exceed[c] + 1) / (N_PERMUTATIONS + 1),
            "occupied_blocks": len(sizes),
            "singleton_blocks": int((sizes == 1).sum()),
        })

    all_results.append({
        "block_size_mb": block_size // 1_000_000,
        "comparison": "RS_gained_vs_all_other_classes",
        "gained_concordance": outcomes[gained].mean(),
        "comparison_concordance": outcomes[~gained].mean(),
        "observed_difference": observed_pooled_difference,
        "permutations": N_PERMUTATIONS,
        "empirical_one_sided_p": (exceed_pooled + 1) / (N_PERMUTATIONS + 1),
        "occupied_blocks": len(sizes),
        "singleton_blocks": int((sizes == 1).sum()),
    })

    print(f"\n{block_size // 1_000_000} Mb blocks:")
    print("  Occupied blocks:", len(sizes))
    print("  Singleton blocks:", int((sizes == 1).sum()))
    print("  Median genes/block:", sizes.median())

results = pd.DataFrame(all_results)

# Correct the four comparisons separately within each block-size analysis.
results["p_adjusted_BH"] = np.nan
for size in results["block_size_mb"].unique():
    mask = results["block_size_mb"] == size
    results.loc[mask, "p_adjusted_BH"] = multipletests(
        results.loc[mask, "empirical_one_sided_p"],
        method="fdr_bh"
    )[1]

out = base + "atac_rnaseq_spatial_block_permutation.csv"
results.to_csv(out, index=False)

print("\n===== SPATIAL-BLOCK PERMUTATION RESULTS =====")
print(results[
    ["block_size_mb", "comparison", "gained_concordance",
     "comparison_concordance", "observed_difference",
     "empirical_one_sided_p", "p_adjusted_BH"]
].to_string(index=False))

print("\nSaved:", out)
print("\nEmpirical p-values use the add-one correction.")
print("These tests preserve class counts within genomic blocks,")
print("but do not eliminate all spatial or biological dependence.")
