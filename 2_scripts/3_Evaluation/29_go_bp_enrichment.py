from pathlib import Path
import gzip
import pandas as pd
from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudy
from statsmodels.stats.multitest import multipletests

ROOT = Path.cwd()
GENE_FILE = ROOT / "3_results/tables/integration/atac_rnaseq_single_class_gene_concordance_v2.csv"
MAP_FILE = ROOT / "1_data_raw/0_reference/hg38/annotation/bed/gene_id_symbol_map.tsv"
GAF_FILE = ROOT / "1_data_raw/4_functional_annotation/GO/HUMAN-uniprot.gaf.gz"
OBO_FILE = ROOT / "1_data_raw/4_functional_annotation/GO/go-basic.obo"
OUT_DIR = ROOT / "3_results/tables/integration/go_enrichment"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_BACKGROUND = 10

# 1. Load single-class genes and uniquely mapped protein-coding symbols
genes = pd.read_csv(GENE_FILE)
mapping = pd.read_csv(MAP_FILE, sep="\t", header=None)
mapping.columns = ["ensembl_id", "ensembl_versioned", "symbol", "gene_type"]

pc = mapping[
    (mapping["gene_type"] == "protein_coding")
    & mapping["symbol"].notna()
].copy()

symbol_counts = pc.groupby("symbol")["ensembl_id"].nunique()
unique_symbols = set(symbol_counts[symbol_counts == 1].index)
pc = pc[pc["symbol"].isin(unique_symbols)]

genes = genes.merge(
    pc[["ensembl_id", "symbol"]],
    on="ensembl_id",
    how="inner",
    validate="one_to_one"
)

# 2. Parse GO ontology
godag = GODag(str(OBO_FILE), optional_attrs={"relationship"})

# 3. Read direct UniProt GO Biological Process annotations.
# Exclude NOT annotations and propagate annotations through GOATOOLS.
symbol_to_go = {}
with gzip.open(GAF_FILE, "rt") as handle:
    for line in handle:
        if line.startswith("!"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 17 or fields[0] != "UniProtKB":
            continue
        if fields[8] != "P" or "NOT" in fields[3].split("|"):
            continue
        symbol = fields[2].strip()
        go_id = fields[4].strip()
        if symbol and go_id in godag and godag[go_id].namespace == "biological_process":
            symbol_to_go.setdefault(symbol, set()).add(go_id)

# Map GO annotations to eligible Ensembl gene IDs.
# Each gene receives direct BP terms; GOATOOLS handles ancestor propagation.
gene_to_go = {}
for row in genes.itertuples(index=False):
    terms = symbol_to_go.get(row.symbol, set())
    if terms:
        gene_to_go[row.ensembl_id] = terms

universe = sorted(gene_to_go)
universe_set = set(universe)

if not universe:
    raise ValueError(
        "No eligible genes mapped to GO BP terms. "
        "Check GAF parsing, symbol mapping, and GO ontology IDs."
    )

# 4. Define gene sets
gained = genes[genes["atac_class"] == "RS_gained"].copy()
concordant = gained[
    gained["concordant"].astype(str).str.lower().isin(["true", "1", "yes"])
].copy()

gene_sets = {
    "RS_gained_all": set(gained["ensembl_id"]) & universe_set,
    "RS_gained_concordant": set(concordant["ensembl_id"]) & universe_set,
}

# 5. Run GOATOOLS with all eligible background genes.
# Use only BP ontology terms; annotations are propagated to ancestors.
go2gene = {}
for gene_id, terms in gene_to_go.items():
    for go_id in terms:
        go2gene.setdefault(go_id, set()).add(gene_id)
        for ancestor in godag[go_id].get_all_parents():
            if ancestor in godag and godag[ancestor].namespace == "biological_process":
                go2gene.setdefault(ancestor, set()).add(gene_id)

# Test every term with >=10 background genes, including zero-overlap terms.
# Apply BH to the full family before display filtering.
rows = []
summary = []

for set_name, selected in gene_sets.items():
    n_selected = len(selected)
    print(f"\n{set_name}: {n_selected} selected genes in GO-annotated background")

    if n_selected == 0:
        print("WARNING: no selected genes in the eligible background.")
        continue

    tested = []
    for go_id, bg_genes in go2gene.items():
        bg_overlap = bg_genes & universe_set
        K = len(bg_overlap)
        if K < MIN_BACKGROUND:
            continue

        overlap = selected & bg_overlap
        k = len(overlap)
        # One-sided hypergeometric over-representation test.
        from scipy.stats import hypergeom
        p_value = float(hypergeom.sf(k - 1, len(universe), K, n_selected))

        tested.append({
            "analysis": set_name,
            "GO_ID": go_id,
            "term": godag[go_id].name,
            "background_n": len(universe),
            "selected_n": n_selected,
            "term_background_n": K,
            "overlap_n": k,
            "expected_overlap": n_selected * K / len(universe),
            "fold_enrichment": (k / n_selected) / (K / len(universe)),
            "p_value": p_value,
            "overlap_ensembl_ids": ";".join(sorted(overlap)),
            "overlap_symbols": ";".join(sorted(
                genes.loc[genes["ensembl_id"].isin(overlap), "symbol"].dropna().unique()
            )),
        })

    if not tested:
        print("WARNING: no GO terms passed the background-size threshold.")
        continue

    adjusted = multipletests(
        [item["p_value"] for item in tested],
        method="fdr_bh"
    )[1]

    for item, q_value in zip(tested, adjusted):
        item["adjusted_p_value"] = float(q_value)
        rows.append(item)

    summary.append({
        "analysis": set_name,
        "analysis_genes_total": len(gained) if set_name == "RS_gained_all" else len(concordant),
        "selected_genes_with_GO_BP": n_selected,
        "background_genes_with_GO_BP": len(universe),
        "GO_terms_tested_BH": len(tested),
        "terms_q_lt_0.05": int((adjusted < 0.05).sum()),
        "terms_q_lt_0.05_and_overlap_ge_3": int(
            sum(q < 0.05 and item["overlap_n"] >= 3
                for q, item in zip(adjusted, tested))
        ),
    })

    print("GO terms tested:", len(tested))
    print("Terms with BH q < 0.05:", int((adjusted < 0.05).sum()))

# 6. Save full tables, plus readable significant-term tables.
result = pd.DataFrame(rows)
if not result.empty:
    result = result.sort_values(
        ["analysis", "adjusted_p_value", "p_value", "overlap_n"],
        ascending=[True, True, True, False]
    )
    result.to_csv(OUT_DIR / "GO_BP_enrichment_all_tested_terms.csv", index=False)

    display = result[
        (result["overlap_n"] >= 3) & (result["adjusted_p_value"] < 0.05)
    ].copy()
    display.to_csv(OUT_DIR / "GO_BP_enrichment_significant_terms.csv", index=False)

    for set_name in gene_sets:
        subset = result[
            (result["analysis"] == set_name) & (result["overlap_n"] >= 3)
        ].head(15)
        print(f"\nTop GO BP terms for {set_name} (displayed after full BH correction):")
        if subset.empty:
            print("No terms with overlap >=3.")
        else:
            print(subset[[
                "GO_ID", "term", "overlap_n", "term_background_n",
                "fold_enrichment", "p_value", "adjusted_p_value"
            ]].to_string(index=False))

pd.DataFrame(summary).to_csv(
    OUT_DIR / "GO_BP_enrichment_run_summary.csv", index=False
)

print("\nCoverage summary:")
print(pd.DataFrame(summary).to_string(index=False))
print("\nSaved outputs to:", OUT_DIR)
