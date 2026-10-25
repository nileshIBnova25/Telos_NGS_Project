#!/usr/bin/env python3
"""R2: Build a descriptive WI-38 young-versus-replicative-senescence RNA table.

The sample/identifier names must first be confirmed in rna_config.json created
by 1_rna_inspect_tpm_table.py. TPM fold changes here are descriptive, not a
replacement for count-based statistical differential-expression testing.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from rna_utils import begin_run, repo_root, resolve_from_repo, setup_logger, write_done_marker


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(root / "1_data_raw" / "4_rnaseq" / "rna_config.json"))
    parser.add_argument("--out-dir", default=str(root / "1_data_raw" / "4_rnaseq" / "outputs"))
    return parser.parse_args()


def fail(logger, message: str) -> int:
    logger.error(message)
    return 2


def main() -> int:
    args = parse_args()
    root = repo_root()
    config_path = resolve_from_repo(args.config, root)
    out_dir = resolve_from_repo(args.out_dir, root)
    base_out_dir = config_path.parent
    log_dir = base_out_dir / "logs"
    done_path = base_out_dir / "R2_rna_groundtruth.done.txt"
    begin_run(done_path)
    logger, log_path = setup_logger("R2_rna_make_groundtruth", log_dir)

    if not config_path.is_file():
        return fail(logger, f"Config not found: {config_path}. Run script 1 first.")
    config = json.loads(config_path.read_text(encoding="utf-8"))

    workbook_value = config.get("workbook_path")
    sheet_name = config.get("sheet_name")
    gene_col = config.get("gene_column")
    young_cols = config.get("young_columns")
    rs_cols = config.get("rs_columns")
    epsilon = float(config.get("epsilon", 0.001))

    if not workbook_value:
        return fail(logger, "Config 'workbook_path' is empty. Rerun script 1 or edit rna_config.json.")
    if not sheet_name:
        return fail(logger, "Config 'sheet_name' is empty.")
    if not gene_col or not isinstance(gene_col, str):
        return fail(logger, "Set 'gene_column' to the exact gene ID/symbol column in rna_config.json.")
    if not isinstance(young_cols, list) or not young_cols:
        return fail(logger, "Set 'young_columns' to the exact WI-38 young/proliferating replicate columns.")
    if not isinstance(rs_cols, list) or not rs_cols:
        return fail(logger, "Set 'rs_columns' to the exact WI-38 PDL46/PDL50 late-passage replicate columns.")
    if not np.isfinite(epsilon) or epsilon <= 0:
        return fail(logger, "Config 'epsilon' must be a positive finite number.")
    if len(set(young_cols)) != len(young_cols) or len(set(rs_cols)) != len(rs_cols):
        return fail(logger, "Duplicate column names were found within a group. Remove duplicates in rna_config.json.")
    overlap = sorted(set(young_cols).intersection(rs_cols))
    if overlap:
        return fail(logger, f"Young and RS groups overlap; remove these columns from one group: {overlap}")

    workbook_path = resolve_from_repo(workbook_value, root)
    if not workbook_path.is_file():
        return fail(logger, f"Workbook not found: {workbook_path}")

    logger.info("Reading workbook=%s sheet=%s", workbook_path, sheet_name)
    df = pd.read_excel(
    workbook_path,
    sheet_name=sheet_name,
    engine="calamine"
)
    logger.info("Loaded dataframe shape=%s", df.shape)

    required = [gene_col, *young_cols, *rs_cols]
    missing = [col for col in required if col not in df.columns]
    if missing:
        logger.error("Columns missing from sheet '%s': %s", sheet_name, missing)
        logger.error("Available columns: %s", df.columns.tolist())
        return 2

    if not df[gene_col].notna().any():
        return fail(logger, f"Gene identifier column '{gene_col}' contains no non-empty values.")

    for col in [*young_cols, *rs_cols]:
        converted = pd.to_numeric(df[col], errors="coerce")
        if converted.notna().sum() == 0:
            return fail(logger, f"Sample column '{col}' has no numeric values; check the selected sheet/column.")
        if (converted.dropna() < 0).any():
            return fail(logger, f"Sample column '{col}' contains negative values, which is unexpected for TPM.")
        n_missing = int(converted.isna().sum())
        if n_missing:
            logger.warning("Column %s has %d blank/non-numeric values; means will skip missing values.", col, n_missing)
        df[col] = converted

    # Means are descriptive summaries across the selected biological replicate columns.
    df["young_tpm_mean"] = df[young_cols].mean(axis=1, skipna=True)
    df["rs_tpm_mean"] = df[rs_cols].mean(axis=1, skipna=True)
    df["young_log2"] = np.log2(df["young_tpm_mean"] + 1.0)
    df["rs_log2"] = np.log2(df["rs_tpm_mean"] + 1.0)
    df["log2FC_rs_vs_young"] = np.log2(
        (df["rs_tpm_mean"] + epsilon) / (df["young_tpm_mean"] + epsilon)
    )

    out = df[[gene_col, "young_tpm_mean", "rs_tpm_mean", "young_log2", "rs_log2", "log2FC_rs_vs_young"]].copy()
    out = out.rename(columns={gene_col: "gene_id_or_symbol"})
    out = out.dropna(subset=["gene_id_or_symbol"])
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "rnaseq_groundtruth_wi38_young_vs_rs.csv"
    out.to_csv(out_path, index=False)

    logger.info("Wrote %s (rows=%d, columns=%d)", out_path, *out.shape)
    logger.info("Young samples (%d): %s", len(young_cols), young_cols)
    logger.info("Late/RS samples (%d): %s", len(rs_cols), rs_cols)
    logger.info("Positive log2FC: %d; negative log2FC: %d", (out.log2FC_rs_vs_young > 0).sum(), (out.log2FC_rs_vs_young < 0).sum())
    write_done_marker(
        done_path,
        "SUCCESS\n"
        f"OUT_CSV={out_path}\nSHEET={sheet_name}\nGENE_COL={gene_col}\n"
        f"YOUNG_COLS={json.dumps(young_cols)}\nRS_COLS={json.dumps(rs_cols)}\nLOG={log_path}",
    )
    logger.info("Wrote completion marker: %s", done_path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
