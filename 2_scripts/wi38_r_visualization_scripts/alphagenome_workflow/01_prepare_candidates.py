#!/usr/bin/env python3
"""Prepare a QC-aware candidate list for WI-38 AlphaGenome locus scans.

Reads the user's ground-truth CSV and RS_TC_raw_TPM worksheet, calculates
per-sample/group expression summaries, and selects up to N genes in each
log2FC direction. It does not treat TPM-derived statistics as definitive DE.
"""
from pathlib import Path
import re
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]  # project root when extracted into 1_data_raw/wi38_r_visualization_scripts/alphagenome_workflow/
# If moved elsewhere, paths below fall back to the current working directory.
PROJECT = ROOT if (ROOT / "1_data_raw").exists() else Path.cwd()
GT_PATH = PROJECT / "1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv"
TPM_PATH = PROJECT / "1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx"
SHEET = "RS_TC_raw_TPM"
OUT_DIR = PROJECT / "1_data_raw/4_rnaseq/outputs/alphagenome"
N_EACH = 20
MIN_TPM_IN_EITHER_GROUP = 1.0

YOUNG_SAMPLES = ["PDL20_TP1_A", "PDL20_TP1_B", "PDL20_TP1_C"]
RS_SAMPLES = ["PDL46_TP6_A", "PDL46_TP6_B", "PDL50_TP7_A", "PDL50_TP7_B"]


def norm(value):
    return re.sub(r"[^a-z0-9]", "", str(value).strip().lower())


def find_col(df, preferred, patterns=(), label="column"):
    normalized = {norm(c): c for c in df.columns}
    for p in preferred:
        if norm(p) in normalized:
            return normalized[norm(p)]
    for c in df.columns:
        n = norm(c)
        if any(re.search(p, n) for p in patterns):
            return c
    raise ValueError(
        f"Cannot detect {label}. Available columns: {list(df.columns)}"
    )


def key(value):
    s = str(value).strip()
    # Remove Ensembl version suffixes but otherwise preserve the identifier.
    s = re.sub(r"\.(\d+)$", "", s)
    return s.upper()


def main():
    print(f"Project root: {PROJECT}")
    print(f"Ground truth: {GT_PATH}")
    print(f"TPM workbook: {TPM_PATH}")
    if not GT_PATH.exists():
        raise FileNotFoundError(f"Ground-truth CSV not found: {GT_PATH}")
    if not TPM_PATH.exists():
        raise FileNotFoundError(f"TPM workbook not found: {TPM_PATH}")

    gt = pd.read_csv(GT_PATH)
    gt.columns = [str(c).strip() for c in gt.columns]
    gt_gene = find_col(
        gt,
        ["gene_id_or_symbol", "gene_symbol", "gene", "gene_name", "symbol", "gene_id"],
        label="gene identifier in ground-truth CSV",
    )
    gt_lfc = find_col(
        gt,
        ["log2FC_rs_vs_young", "log2FC", "log2_fold_change", "log2FoldChange"],
        patterns=(r"log2.*fc", r"log2.*foldchange", r"foldchange.*log2"),
        label="log2 fold-change column",
    )
    print(f"Ground-truth gene column: {gt_gene!r}; effect column: {gt_lfc!r}")

    # python-calamine supports Strict Open XML files that may appear as zero-sheet
    # workbooks to openpyxl. Install it in the active environment if needed.
    tpm = pd.read_excel(TPM_PATH, sheet_name=SHEET, engine="calamine")
    tpm.columns = [str(c).strip() for c in tpm.columns]
    tpm_gene = find_col(
        tpm,
        ["gene", "gene_symbol", "gene_name", "symbol", "gene_id", "gene_id_or_symbol"],
        patterns=(r"^gene$", r"genesymbol", r"genename", r"geneid"),
        label="gene identifier in TPM sheet",
    )
    cols_norm = {norm(c): c for c in tpm.columns}
    missing = [s for s in YOUNG_SAMPLES + RS_SAMPLES if norm(s) not in cols_norm]
    if missing:
        raise ValueError(
            "Expected sample columns were not found: " + str(missing) +
            "\nActual TPM columns: " + str(list(tpm.columns))
        )
    young_cols = [cols_norm[norm(s)] for s in YOUNG_SAMPLES]
    rs_cols = [cols_norm[norm(s)] for s in RS_SAMPLES]

    gt = gt[[gt_gene, gt_lfc] + [c for c in gt.columns if c not in (gt_gene, gt_lfc)]].copy()
    gt["_key"] = gt[gt_gene].map(key)
    gt["ground_truth_log2FC_RS_vs_Young"] = pd.to_numeric(gt[gt_lfc], errors="coerce")
    gt = gt.dropna(subset=["_key", "ground_truth_log2FC_RS_vs_Young"])
    gt = gt.drop_duplicates("_key", keep="first")

    keep_cols = [tpm_gene] + young_cols + rs_cols
    exp = tpm[keep_cols].copy()
    exp["_key"] = exp[tpm_gene].map(key)
    for c in young_cols + rs_cols:
        exp[c] = pd.to_numeric(exp[c], errors="coerce")
    exp = exp.groupby("_key", as_index=False).agg(
        {tpm_gene: "first", **{c: "mean" for c in young_cols + rs_cols}}
    )
    merged = gt.merge(exp, on="_key", how="left", suffixes=("", "_tpm"))
    match_rate = merged[young_cols + rs_cols].notna().any(axis=1).mean()
    print(f"Ground-truth rows: {len(gt):,}; TPM key match rate: {match_rate:.1%}")
    if match_rate < 0.10:
        print("Example ground-truth identifiers:", gt[gt_gene].head(10).tolist())
        print("Example TPM identifiers:", tpm[tpm_gene].head(10).tolist())
        raise ValueError(
            "Fewer than 10% of ground-truth identifiers matched the TPM gene column. "
            "The two files may use different ID types; inspect the examples above before proceeding."
        )

    merged["mean_TPM_Young"] = merged[young_cols].mean(axis=1, skipna=True)
    merged["mean_TPM_RS"] = merged[rs_cols].mean(axis=1, skipna=True)
    merged["max_group_mean_TPM"] = merged[["mean_TPM_Young", "mean_TPM_RS"]].max(axis=1)
    merged["TPM_log2_mean_ratio_RS_vs_Young"] = np.log2(
        (merged["mean_TPM_RS"].fillna(0) + 0.1) /
        (merged["mean_TPM_Young"].fillna(0) + 0.1)
    )
    merged["effect_direction_agrees_with_group_means"] = (
        np.sign(merged["ground_truth_log2FC_RS_vs_Young"]) ==
        np.sign(merged["TPM_log2_mean_ratio_RS_vs_Young"])
    )
    young_median = merged[young_cols].median(axis=1, skipna=True)
    rs_median = merged[rs_cols].median(axis=1, skipna=True)
    merged["n_RS_above_Young_median"] = merged[rs_cols].gt(young_median, axis=0).sum(axis=1)
    merged["n_RS_below_Young_median"] = merged[rs_cols].lt(young_median, axis=0).sum(axis=1)
    merged["n_Young_below_RS_median"] = merged[young_cols].lt(rs_median, axis=0).sum(axis=1)
    merged["n_Young_above_RS_median"] = merged[young_cols].gt(rs_median, axis=0).sum(axis=1)
    merged["direction"] = np.where(
        merged["ground_truth_log2FC_RS_vs_Young"] > 0, "Up in RS", "Down in RS"
    )
    merged["abs_ground_truth_log2FC"] = merged["ground_truth_log2FC_RS_vs_Young"].abs()

    # Only use genes with adequate TPM representation for this first sequence scan;
    # don't discard direction disagreements—flag and rank them lower instead.
    eligible = merged[
        (merged["max_group_mean_TPM"] >= MIN_TPM_IN_EITHER_GROUP) &
        (merged["ground_truth_log2FC_RS_vs_Young"] != 0)
    ].copy()
    eligible["consistency_score"] = np.where(
        eligible["direction"] == "Up in RS",
        (eligible["n_RS_above_Young_median"] / len(rs_cols) +
         eligible["n_Young_below_RS_median"] / len(young_cols)) / 2,
        (eligible["n_RS_below_Young_median"] / len(rs_cols) +
         eligible["n_Young_above_RS_median"] / len(young_cols)) / 2,
    )
    eligible = eligible.sort_values(
        ["direction", "effect_direction_agrees_with_group_means", "consistency_score", "abs_ground_truth_log2FC", "max_group_mean_TPM"],
        ascending=[True, False, False, False, False],
    )
    candidate_parts = []
    for direction in ["Up in RS", "Down in RS"]:
        part = eligible[eligible["direction"] == direction].head(N_EACH).copy()
        candidate_parts.append(part)
        print(f"{direction}: selected {len(part)} / requested {N_EACH}")
    candidates = pd.concat(candidate_parts, ignore_index=True)
    # Standardize identifier columns. When both source tables use the same column
    # name pandas suffixes the TPM identifier with "_tpm" during the merge.
    merged_tpm_gene = tpm_gene
    if merged_tpm_gene == gt_gene and f"{tpm_gene}_tpm" in candidates.columns:
        merged_tpm_gene = f"{tpm_gene}_tpm"
    rename_map = {}
    if merged_tpm_gene in candidates.columns:
        rename_map[merged_tpm_gene] = "tpm_gene_identifier"
    if gt_gene in candidates.columns:
        rename_map[gt_gene] = "ground_truth_gene_identifier"
    candidates.rename(columns=rename_map, inplace=True)
    if "tpm_gene_identifier" not in candidates.columns:
        candidates["tpm_gene_identifier"] = candidates["_key"]
    candidates["gene_for_lookup"] = candidates["tpm_gene_identifier"].astype(str)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_path = OUT_DIR / "candidate_gene_scores_all.csv"
    out_path = OUT_DIR / "candidate_genes_top20_each_direction.csv"
    merged.to_csv(all_path, index=False)
    candidates.to_csv(out_path, index=False)
    print(f"\nWrote all joined gene summaries: {all_path}")
    print(f"Wrote AlphaGenome shortlist: {out_path}")
    print("\nSelected candidates:")
    show = ["gene_for_lookup", "direction", "ground_truth_log2FC_RS_vs_Young", "mean_TPM_Young", "mean_TPM_RS", "effect_direction_agrees_with_group_means", "consistency_score"]
    show = [c for c in show if c in candidates.columns]
    print(candidates[show].to_string(index=False))
    print("\nReminder: TPM-based consistency summaries are QC aids, not count-based differential-expression tests.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
