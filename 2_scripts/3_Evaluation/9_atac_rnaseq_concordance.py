import os
import pandas as pd
import numpy as np

candidate_file = "3_results/tables/integrated_candidate_regions.csv"
integration_file = (
    "3_results/tables/integration/"
    "atac_region_gene_rnaseq_associations.csv"
)
out_dir = "3_results/tables/integration/"
os.makedirs(out_dir, exist_ok=True)

candidates = pd.read_csv(candidate_file)
integrated = pd.read_csv(integration_file)

# Validate required columns.
required_candidates = {
    "name", "label", "log2FC_rs_vs_young"
}
required_integrated = {
    "region_name", "ensembl_id", "gene_symbol",
    "feature_type", "rna_match_status",
    "log2FoldChange", "padj"
}

if not required_candidates.issubset(candidates.columns):
    raise ValueError(
        f"Missing candidate columns: "
        f"{required_candidates - set(candidates.columns)}"
    )

if not required_integrated.issubset(integrated.columns):
    raise ValueError(
        f"Missing integration columns: "
        f"{required_integrated - set(integrated.columns)}"
    )

if candidates["name"].duplicated().any():
    raise ValueError("Candidate region names are not unique.")

# Check whether ATAC class labels agree with the observed ATAC effect.
candidates["observed_atac_direction"] = np.select(
    [
        candidates["log2FC_rs_vs_young"] > 0,
        candidates["log2FC_rs_vs_young"] < 0
    ],
    ["positive", "negative"],
    default="zero"
)

print("===== ATAC class versus observed effect =====")
print(
    pd.crosstab(
        candidates["label"],
        candidates["observed_atac_direction"],
        margins=True
    ).to_string()
)

# Join RNA associations to the observed ATAC effect by region ID.
region_info = candidates[
    ["name", "label", "log2FC_rs_vs_young"]
].rename(columns={
    "name": "region_name",
    "label": "atac_class",
    "log2FC_rs_vs_young": "observed_atac_log2FC"
})

df = integrated.drop(columns=["atac_class"], errors="ignore").merge(
    region_info,
    on="region_name",
    how="left",
    validate="many_to_one"
)

if df["observed_atac_log2FC"].isna().any():
    missing = df.loc[
        df["observed_atac_log2FC"].isna(), "region_name"
    ].nunique()
    raise ValueError(
        f"{missing} integrated regions lack candidate-table matches."
    )

# Primary analysis: unambiguous symbol matches with significant RNA results.
eligible = df[
    (df["rna_match_status"] == "matched_unique_symbol")
    & df["padj"].notna()
    & (df["padj"] < 0.05)
    & df["log2FoldChange"].notna()
].copy()

# Avoid counting the same region-gene pair twice because it overlaps
# both a gene body and a promoter.
eligible = eligible.drop_duplicates(
    subset=["region_name", "ensembl_id"]
).copy()

eligible["atac_direction"] = np.select(
    [
        eligible["observed_atac_log2FC"] > 0,
        eligible["observed_atac_log2FC"] < 0
    ],
    ["positive", "negative"],
    default="zero"
)

eligible["rna_direction"] = np.select(
    [
        eligible["log2FoldChange"] > 0,
        eligible["log2FoldChange"] < 0
    ],
    ["positive", "negative"],
    default="zero"
)

eligible["concordance"] = np.select(
    [
        (eligible["atac_direction"] == "zero")
        | (eligible["rna_direction"] == "zero"),
        eligible["atac_direction"] == eligible["rna_direction"]
    ],
    ["zero_effect", "concordant"],
    default="discordant"
)

# Save association-level results.
association_path = (
    out_dir + "atac_rnaseq_direction_concordance_associations.csv"
)
eligible.to_csv(association_path, index=False)

# Summarize by class; each region-gene pair is counted once.
summary = eligible.groupby(
    ["atac_class", "concordance"], dropna=False
).agg(
    region_gene_pairs=("region_name", "size"),
    distinct_regions=("region_name", "nunique"),
    distinct_genes=("ensembl_id", "nunique")
).reset_index()

summary_path = (
    out_dir + "atac_rnaseq_direction_concordance_summary.csv"
)
summary.to_csv(summary_path, index=False)

# Calculate percentages among non-zero, directionally classifiable pairs.
valid = eligible[eligible["concordance"] != "zero_effect"].copy()

if len(valid):
    class_percentages = valid.groupby("atac_class").agg(
        classifiable_pairs=("concordance", "size"),
        concordant_pairs=(
            "concordance",
            lambda x: (x == "concordant").sum()
        ),
        discordant_pairs=(
            "concordance",
            lambda x: (x == "discordant").sum()
        )
    ).reset_index()

    class_percentages["concordance_percent"] = (
        100 * class_percentages["concordant_pairs"]
        / class_percentages["classifiable_pairs"]
    )
else:
    class_percentages = pd.DataFrame()

percentage_path = (
    out_dir + "atac_rnaseq_direction_concordance_percentages.csv"
)
class_percentages.to_csv(percentage_path, index=False)

# Gene-level summary: average observed ATAC effect across distinct
# regions linked to each Ensembl gene, then compare with its RNA effect.
gene_level = eligible.groupby(
    ["ensembl_id", "gene_symbol"], dropna=False
).agg(
    atac_class=("atac_class", lambda x: x.mode().iloc[0]),
    distinct_regions=("region_name", "nunique"),
    mean_observed_atac_log2FC=("observed_atac_log2FC", "mean"),
    rna_log2FoldChange=("log2FoldChange", "first"),
    rna_padj=("padj", "first")
).reset_index()

gene_level["atac_direction"] = np.select(
    [
        gene_level["mean_observed_atac_log2FC"] > 0,
        gene_level["mean_observed_atac_log2FC"] < 0
    ],
    ["positive", "negative"],
    default="zero"
)

gene_level["rna_direction"] = np.select(
    [
        gene_level["rna_log2FoldChange"] > 0,
        gene_level["rna_log2FoldChange"] < 0
    ],
    ["positive", "negative"],
    default="zero"
)

gene_level["concordance"] = np.select(
    [
        (gene_level["atac_direction"] == "zero")
        | (gene_level["rna_direction"] == "zero"),
        gene_level["atac_direction"] == gene_level["rna_direction"]
    ],
    ["zero_effect", "concordant"],
    default="discordant"
)

gene_path = out_dir + "atac_rnaseq_direction_concordance_genes.csv"
gene_level.to_csv(gene_path, index=False)

print("\n===== Primary analysis =====")
print("Significant, uniquely matched region-gene pairs:", len(eligible))
print("Distinct regions:", eligible["region_name"].nunique())
print("Distinct Ensembl genes:", eligible["ensembl_id"].nunique())
print("Concordant pairs:", (eligible["concordance"] == "concordant").sum())
print("Discordant pairs:", (eligible["concordance"] == "discordant").sum())
print("Pairs with a zero effect:", (eligible["concordance"] == "zero_effect").sum())

print("\n===== Concordance summary by ATAC class =====")
print(summary.to_string(index=False))

print("\n===== Concordance percentages =====")
if len(class_percentages):
    print(class_percentages.to_string(index=False))
else:
    print("No pairs available for percentage calculation.")

print("\n===== Gene-level summary =====")
print(
    gene_level.groupby(["atac_class", "concordance"])
    .agg(distinct_genes=("ensembl_id", "nunique"))
    .to_string()
)

print("\nSaved:")
print(association_path)
print(summary_path)
print(percentage_path)
print(gene_path)
