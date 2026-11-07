import os
import pandas as pd

base = "3_results/tables/atac_annotation/"
mapping_file = (
    "1_data_raw/0_reference/hg38/annotation/bed/"
    "gene_id_symbol_map.tsv"
)
rna_file = (
    "3_results/tables/rnaseq/"
    "RS_DESeq2_timecourse_batch_adjusted_unshrunken.csv"
)
candidate_file = "3_results/tables/integrated_candidate_regions.csv"
out_dir = "3_results/tables/integration/"

os.makedirs(out_dir, exist_ok=True)

# Load GENCODE mapping, DESeq2 results, and candidate-region metadata.
mapping = pd.read_csv(
    mapping_file,
    sep="\t",
    header=None,
    names=[
        "ensembl_id",
        "ensembl_id_versioned",
        "gene_symbol",
        "gene_type_mapping"
    ],
    dtype=str
)

rna = pd.read_csv(rna_file)
candidates = pd.read_csv(candidate_file)

# Check that the expected identifiers are available.
required_rna = {"gene", "log2FoldChange", "padj"}
if not required_rna.issubset(rna.columns):
    raise ValueError(
        "RNA-seq file is missing columns: "
        + str(required_rna - set(rna.columns))
    )

if mapping["ensembl_id"].duplicated().any():
    raise ValueError("Duplicate Ensembl IDs found in the GENCODE mapping.")

# Count Ensembl IDs associated with each gene symbol.
# Do not remove ambiguous mappings: preserve them and flag them.
symbol_counts = (
    mapping.dropna(subset=["gene_symbol"])
    .groupby("gene_symbol")["ensembl_id"]
    .nunique()
)

mapping["symbol_mapping_count"] = (
    mapping["gene_symbol"].map(symbol_counts).fillna(0).astype(int)
)

mapping["mapping_status"] = "missing_symbol"
mapping.loc[
    mapping["symbol_mapping_count"] == 1, "mapping_status"
] = "unique"
mapping.loc[
    mapping["symbol_mapping_count"] > 1, "mapping_status"
] = "ambiguous"

# Prepare BEDTools intersection files.
# Each row has 5 benchmark columns followed by 7 annotation columns.
def read_intersections(path, feature):
    df = pd.read_csv(path, sep="\t", header=None)

    df.columns = [
        "chr", "start", "end", "region_name", "atac_class",
        "annotation_chr", "annotation_start", "annotation_end",
        "ensembl_id_versioned", "annotation_gene_symbol",
        "strand", "gene_type_annotation"
    ]

    df["feature_type"] = feature
    df["ensembl_id"] = (
        df["ensembl_id_versioned"]
        .astype(str)
        .str.replace(r"\.\d+$", "", regex=True)
    )
    return df


body = read_intersections(
    base + "atac_gene_body_overlaps.tsv", "gene_body"
)
promoter = read_intersections(
    base + "atac_promoter_overlaps.tsv", "promoter"
)

overlaps = pd.concat([body, promoter], ignore_index=True)

# Attach the complete GENCODE mapping by Ensembl ID.
overlaps = overlaps.merge(
    mapping[
        [
            "ensembl_id",
            "gene_symbol",
            "gene_type_mapping",
            "symbol_mapping_count",
            "mapping_status"
        ]
    ],
    on="ensembl_id",
    how="left",
    validate="many_to_one"
)

# Explicitly identify annotation records without a mapping.
overlaps["mapping_status"] = overlaps["mapping_status"].fillna(
    "ensembl_id_not_in_mapping"
)

# Use RNA-seq gene symbols only where the symbol is unambiguous
# in GENCODE. This avoids assigning one symbol-level result to
# multiple Ensembl genes.
rna_symbols = set(rna["gene"].dropna().astype(str))
overlaps["gene_symbol_match"] = (
    overlaps["gene_symbol"].notna()
    & overlaps["gene_symbol"].isin(rna_symbols)
)

overlaps["rna_match_status"] = "missing_annotation_symbol"

overlaps.loc[
    overlaps["mapping_status"] == "ensembl_id_not_in_mapping",
    "rna_match_status"
] = "ensembl_id_not_in_mapping"

overlaps.loc[
    overlaps["mapping_status"] == "ambiguous",
    "rna_match_status"
] = "ambiguous_symbol_mapping"

unique_symbol = overlaps["mapping_status"] == "unique"

overlaps.loc[
    unique_symbol & overlaps["gene_symbol_match"],
    "rna_match_status"
] = "matched_unique_symbol"

overlaps.loc[
    unique_symbol & ~overlaps["gene_symbol_match"],
    "rna_match_status"
] = "symbol_not_in_rna_results"

# Exclude ambiguous symbols from the DESeq2 join.
# Keep their annotation rows in the output, with RNA statistics NA.
rna_unique = rna.dropna(subset=["gene"]).copy()
rna_unique["gene"] = rna_unique["gene"].astype(str)

if rna_unique["gene"].duplicated().any():
    duplicates = rna_unique.loc[
        rna_unique["gene"].duplicated(keep=False), "gene"
    ].unique()
    raise ValueError(
        "Duplicate gene symbols in DESeq2 results: "
        + ", ".join(duplicates[:10])
    )

rna_unique = rna_unique.rename(columns={"gene": "rna_gene_symbol"})

overlaps = overlaps.merge(
    rna_unique,
    left_on=[
        overlaps["gene_symbol"].where(unique_symbol)
    ],
    right_on="rna_gene_symbol",
    how="left",
    validate="many_to_one"
)

# The join above uses a Series as the left key; standardize the
# resulting symbol-match field and retain the full annotation symbol.
if "rna_gene_symbol" not in overlaps.columns:
    overlaps["rna_gene_symbol"] = pd.NA

# Significance is evaluated only where an RNA-seq result was matched.
overlaps["rna_significant_padj_0.05"] = (
    overlaps["rna_match_status"].eq("matched_unique_symbol")
    & overlaps["padj"].notna()
    & (overlaps["padj"] < 0.05)
)

overlaps["rna_direction"] = "unavailable_or_not_significant"
overlaps.loc[
    overlaps["rna_significant_padj_0.05"]
    & (overlaps["log2FoldChange"] > 0),
    "rna_direction"
] = "positive_timecourse_effect"

overlaps.loc[
    overlaps["rna_significant_padj_0.05"]
    & (overlaps["log2FoldChange"] < 0),
    "rna_direction"
] = "negative_timecourse_effect"

# Save all region-gene-feature associations.
full_path = out_dir + "atac_region_gene_rnaseq_associations.csv"
overlaps.to_csv(full_path, index=False)

# Summary by ATAC class and genomic feature.
summary = overlaps.groupby(
    ["atac_class", "feature_type"], dropna=False
).agg(
    overlap_pairs=("region_name", "size"),
    unique_regions=("region_name", "nunique"),
    unique_ensembl_genes=("ensembl_id", "nunique"),
    unique_symbols=("gene_symbol", "nunique"),
    ambiguous_mapping_rows=(
        "mapping_status",
        lambda x: (x == "ambiguous").sum()
    ),
    unique_symbol_rna_matches=(
        "rna_match_status",
        lambda x: (x == "matched_unique_symbol").sum()
    ),
    significant_rna_association_rows=(
        "rna_significant_padj_0.05", "sum"
    )
).reset_index()

summary_path = out_dir + "atac_rnaseq_integration_summary.csv"
summary.to_csv(summary_path, index=False)

# Separate mapping-status audit.
mapping_summary = overlaps.groupby(
    ["mapping_status", "rna_match_status"], dropna=False
).agg(
    association_rows=("region_name", "size"),
    unique_regions=("region_name", "nunique"),
    unique_ensembl_ids=("ensembl_id", "nunique")
).reset_index()

mapping_summary_path = (
    out_dir + "atac_rnaseq_mapping_status_summary.csv"
)
mapping_summary.to_csv(mapping_summary_path, index=False)

print("===== Integration validation =====")
print("Candidate regions:", candidates["name"].nunique())
print("Region-gene-feature association rows:", len(overlaps))
print("Unique ATAC regions represented:", overlaps["region_name"].nunique())
print("Unique Ensembl IDs represented:", overlaps["ensembl_id"].nunique())
print(
    "Rows with a unique-symbol RNA-seq match:",
    (overlaps["rna_match_status"] == "matched_unique_symbol").sum()
)
print(
    "Rows with ambiguous symbol mappings:",
    (overlaps["rna_match_status"] == "ambiguous_symbol_mapping").sum()
)
print(
    "Rows with symbols absent from RNA-seq:",
    (overlaps["rna_match_status"] == "symbol_not_in_rna_results").sum()
)
print(
    "Significant matched RNA-seq association rows:",
    overlaps["rna_significant_padj_0.05"].sum()
)

print("\n===== RNA-seq mapping status =====")
print(mapping_summary.to_string(index=False))

print("\n===== Summary by ATAC class and feature =====")
print(summary.to_string(index=False))

print("\nSaved:")
print(full_path)
print(summary_path)
print(mapping_summary_path)
