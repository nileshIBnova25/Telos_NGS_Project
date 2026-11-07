import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

input_file = (
    "3_results/tables/integration/"
    "atac_rnaseq_direction_concordance_associations.csv"
)
out_dir = "3_results/tables/integration/"
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(input_file)

required = [
    "region_name", "ensembl_id", "gene_symbol", "atac_class",
    "observed_atac_log2FC", "log2FoldChange", "rna_match_status",
    "padj"
]
missing = [col for col in required if col not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Retain significant RNA results with unambiguous gene mapping.
df = df[
    (df["rna_match_status"] == "matched_unique_symbol")
    & (df["padj"] < 0.05)
    & df["log2FoldChange"].notna()
    & df["observed_atac_log2FC"].notna()
    & df["ensembl_id"].notna()
].copy()

# One row per region-gene pair.
df = df.drop_duplicates(["region_name", "ensembl_id"])

# Keep genes linked to exactly one ATAC class.
gene_classes = df.groupby("ensembl_id")["atac_class"].nunique()
single_class_ids = gene_classes[gene_classes == 1].index
df = df[df["ensembl_id"].isin(single_class_ids)].copy()

# Aggregate multiple regions linked to the same gene within its class.
genes = df.groupby(
    ["ensembl_id", "gene_symbol", "atac_class"], as_index=False
).agg(
    n_regions=("region_name", "nunique"),
    mean_atac_log2FC=("observed_atac_log2FC", "mean"),
    rna_log2FC=("log2FoldChange", "first")
)

genes["concordant"] = (
    np.sign(genes["mean_atac_log2FC"])
    == np.sign(genes["rna_log2FC"])
)
genes = genes[
    (genes["mean_atac_log2FC"] != 0)
    & (genes["rna_log2FC"] != 0)
].copy()

association_out = os.path.join(
    out_dir, "atac_rnaseq_single_class_gene_concordance_v2.csv"
)
genes.to_csv(association_out, index=False)

# Per-class rate and Wilson 95% confidence interval.
rows = []
for cls, group in genes.groupby("atac_class"):
    n = len(group)
    k = int(group["concordant"].sum())
    low, high = proportion_confint(k, n, method="wilson")
    rows.append({
        "atac_class": cls,
        "genes": n,
        "concordant": k,
        "discordant": n - k,
        "concordance_percent": 100 * k / n if n else np.nan,
        "ci95_lower_percent": 100 * low,
        "ci95_upper_percent": 100 * high
    })

summary = pd.DataFrame(rows).sort_values("atac_class")
summary_out = os.path.join(
    out_dir, "atac_rnaseq_single_class_gene_concordance_v2_summary.csv"
)
summary.to_csv(summary_out, index=False)

# Compare RS_gained with each other class using a 2x2 Fisher exact test.
tests = []
if "RS_gained" in set(genes["atac_class"]):
    gained = genes[genes["atac_class"] == "RS_gained"]
    for cls in sorted(set(genes["atac_class"]) - {"RS_gained"}):
        other = genes[genes["atac_class"] == cls]
        table = [
            [int(gained["concordant"].sum()),
             int((~gained["concordant"]).sum())],
            [int(other["concordant"].sum()),
             int((~other["concordant"]).sum())]
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        tests.append({
            "comparison": f"RS_gained_vs_{cls}",
            "gained_genes": len(gained),
            "gained_concordant": int(gained["concordant"].sum()),
            "other_genes": len(other),
            "other_concordant": int(other["concordant"].sum()),
            "odds_ratio": odds_ratio,
            "p_value": p_value
        })

tests = pd.DataFrame(tests)
if not tests.empty:
    tests["p_adjusted_BH"] = multipletests(
        tests["p_value"], method="fdr_bh"
    )[1]

tests_out = os.path.join(
    out_dir, "atac_rnaseq_single_class_gene_concordance_v2_tests.csv"
)
tests.to_csv(tests_out, index=False)

print("===== SINGLE-CLASS GENE CONCORDANCE =====")
print("Genes linked to only one ATAC class before direction filtering:",
      len(single_class_ids))
print("Directionally classifiable genes:", len(genes))
print(summary.to_string(index=False))

print("\n===== RS_GAINED VERSUS OTHER CLASSES =====")
if tests.empty:
    print("No comparisons available.")
else:
    print(tests.to_string(index=False))

print("\nSaved:")
print(association_out)
print(summary_out)
print(tests_out)
print("\nInterpretation: exploratory gene-level sensitivity analysis.")
print("Fisher tests do not account for all biological dependence between genes.")
