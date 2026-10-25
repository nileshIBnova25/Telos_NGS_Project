#!/usr/bin/env python3
"""R3: Create RNA QC summary and exploratory top up/down gene tables."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

from rna_utils import begin_run, repo_root, resolve_from_repo, setup_logger, write_done_marker

REQUIRED = ["gene_id_or_symbol", "young_tpm_mean", "rs_tpm_mean", "young_log2", "rs_log2", "log2FC_rs_vs_young"]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", default=str(root / "1_data_raw" / "4_rnaseq" / "outputs" / "rnaseq_groundtruth_wi38_young_vs_rs.csv")
    )
    parser.add_argument("--out-dir", default=str(root / "1_data_raw" / "4_rnaseq" / "qc"))
    parser.add_argument("--top-n", type=int, default=50, help="Number of genes in each exploratory top list (default: 50).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = repo_root()
    input_path = resolve_from_repo(args.input, root)
    out_dir = resolve_from_repo(args.out_dir, root)
    base_out_dir = root / "1_data_raw" / "4_rnaseq"
    log_dir = base_out_dir / "logs"
    done_path = base_out_dir / "R3_rna_qc.done.txt"
    begin_run(done_path)
    logger, log_path = setup_logger("R3_rna_qc_summary", log_dir)

    if args.top_n < 1:
        logger.error("--top-n must be >= 1")
        return 2
    if not input_path.is_file():
        logger.error("Ground-truth CSV not found: %s", input_path)
        logger.error("Run 2_rna_make_groundtruth.py successfully first.")
        return 2

    df = pd.read_csv(input_path)
    logger.info("Loaded %s shape=%s", input_path, df.shape)
    missing = [column for column in REQUIRED if column not in df.columns]
    if missing:
        logger.error("Input is missing required columns: %s", missing)
        logger.error("Available columns: %s", df.columns.tolist())
        return 2

    numeric_columns = REQUIRED[1:]
    for col in numeric_columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    n_total = int(len(df))
    qc = {
        "input_file": str(input_path),
        "n_genes_total": n_total,
        "n_genes_young_tpm_gt1": int((df["young_tpm_mean"] > 1).sum()),
        "n_genes_rs_tpm_gt1": int((df["rs_tpm_mean"] > 1).sum()),
        "n_genes_tpm_gt1_both_groups": int(((df["young_tpm_mean"] > 1) & (df["rs_tpm_mean"] > 1)).sum()),
        "n_genes_zero_tpm_both_groups": int(((df["young_tpm_mean"] == 0) & (df["rs_tpm_mean"] == 0)).sum()),
        "n_missing_log2FC": int(df["log2FC_rs_vs_young"].isna().sum()),
        "median_young_tpm": float(df["young_tpm_mean"].median()),
        "median_rs_tpm": float(df["rs_tpm_mean"].median()),
        "median_log2FC": float(df["log2FC_rs_vs_young"].median()),
        "n_log2FC_positive": int((df["log2FC_rs_vs_young"] > 0).sum()),
        "n_log2FC_negative": int((df["log2FC_rs_vs_young"] < 0).sum()),
        "n_log2FC_zero": int((df["log2FC_rs_vs_young"] == 0).sum()),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_path = out_dir / "rnaseq_qc_summary.csv"
    pd.DataFrame([qc]).to_csv(qc_path, index=False)

    valid = df.dropna(subset=["log2FC_rs_vs_young"]).copy()
    top_up = valid.sort_values("log2FC_rs_vs_young", ascending=False, kind="mergesort").head(args.top_n)
    top_down = valid.sort_values("log2FC_rs_vs_young", ascending=True, kind="mergesort").head(args.top_n)
    top_up_path = out_dir / "top_up_genes.csv"
    top_down_path = out_dir / "top_down_genes.csv"
    top_up.to_csv(top_up_path, index=False)
    top_down.to_csv(top_down_path, index=False)

    logger.info("Wrote QC summary: %s", qc_path)
    logger.info("Wrote top up genes: %s (n=%d)", top_up_path, len(top_up))
    logger.info("Wrote top down genes: %s (n=%d)", top_down_path, len(top_down))
    logger.warning("Top-gene lists are ranked by descriptive TPM log2FC; they are not significance-filtered DE results.")
    write_done_marker(
        done_path,
        f"SUCCESS\nQC_SUMMARY={qc_path}\nTOP_UP={top_up_path}\nTOP_DOWN={top_down_path}\nLOG={log_path}",
    )
    logger.info("Wrote completion marker: %s", done_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
