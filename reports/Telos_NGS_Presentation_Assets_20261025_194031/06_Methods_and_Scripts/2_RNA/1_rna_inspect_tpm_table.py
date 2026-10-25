
#!/usr/bin/env python3
"""
R1: Inspect the GEO GSE175533 TPM workbook and prepare an RNA config template.

Run from anywhere:
    python 2_scripts/2_RNA/1_rna_inspect_tpm_table.py

Uses pandas with the Calamine Excel engine to inspect worksheet names and
headers. Does not guess sample groups.

Outputs are stored locally under:
    1_data_raw/4_rnaseq/inspection/
    1_data_raw/4_rnaseq/rna_config.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from rna_utils import (
    begin_run,
    repo_root,
    resolve_from_repo,
    setup_logger,
    write_done_marker,
)


DEFAULT_FILENAME = "GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx"
DEFAULT_SHEET = "RS_TC_raw_TPM"

KEYWORDS = [
    "PDL20",
    "PDL46",
    "PDL50",
    "TP1",
    "TP6",
    "TP7",
    "WI-38",
    "WI38",
    "RS",
    "hTERT",
    "RIS",
    "CD",
    "young",
    "senesc",
    "replicate",
]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    root = repo_root()

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "--input",
        default=str(root / "1_data_raw" / DEFAULT_FILENAME),
        help=(
            "Path to the TPM Excel workbook. "
            "Defaults to 1_data_raw/<workbook>."
        ),
    )

    parser.add_argument(
        "--sheet",
        default=DEFAULT_SHEET,
        help="Exact name of the main RS time-course worksheet.",
    )

    parser.add_argument(
        "--out-dir",
        default=str(root / "1_data_raw" / "4_rnaseq"),
        help="Directory for local inspection reports and RNA configuration.",
    )

    return parser.parse_args()


def read_workbook_headers(
    input_path: Path,
    logger,
) -> tuple[list[dict], list[dict], dict[str, list[str]]]:
    """
    Read worksheet names and column headers using Calamine.

    Only the header row is read from each worksheet; the full expression
    matrices are not loaded into memory.
    """
    sheet_rows: list[dict] = []
    column_rows: list[dict] = []
    all_headers: dict[str, list[str]] = {}

    with pd.ExcelFile(input_path, engine="calamine") as workbook:
        sheet_names = workbook.sheet_names

        logger.info("Workbook sheets: %s", sheet_names)

        if not sheet_names:
            raise ValueError(
                "The Excel reader found no worksheets. "
                "Check whether the workbook is valid and readable."
            )

        for sheet_name in sheet_names:
            header_df = pd.read_excel(
                workbook,
                sheet_name=sheet_name,
                header=0,
                nrows=0,
            )

            headers = [
                str(value).strip() if value is not None else ""
                for value in header_df.columns.tolist()
            ]

            all_headers[sheet_name] = headers

            sheet_rows.append(
                {
                    "sheet_name": sheet_name,
                    # Reading only headers means total data rows are unknown.
                    "reported_max_row": "",
                    "reported_max_column": len(headers),
                    "header_count": len(headers),
                }
            )

            for position, column_name in enumerate(headers, start=1):
                column_rows.append(
                    {
                        "sheet_name": sheet_name,
                        "column_position_1based": position,
                        "column_name": column_name,
                    }
                )

            logger.info(
                "Inspected sheet '%s': %d columns",
                sheet_name,
                len(headers),
            )

    return sheet_rows, column_rows, all_headers


def write_config_template(
    config_path: Path,
    input_path: Path,
    root: Path,
    sheet_name: str,
    logger,
) -> None:
    """Create a local RNA configuration template without overwriting it."""
    if config_path.exists():
        logger.warning(
            "Config already exists; leaving it unchanged: %s",
            config_path,
        )
        return

    try:
        workbook_path = str(input_path.relative_to(root))
    except ValueError:
        workbook_path = str(input_path)

    config = {
        "workbook_path": workbook_path,
        "sheet_name": sheet_name,
        "gene_column": None,
        "young_columns": [],
        "rs_columns": [],
        "epsilon": 0.001,
        "notes": [
            (
                "Fill in gene_column with the exact gene identifier column "
                "from the selected worksheet."
            ),
            (
                "Fill young_columns with the confirmed WI-38 PDL20 "
                "replicate column names."
            ),
            (
                "Fill rs_columns with the confirmed WI-38 PDL46 and PDL50 "
                "replicate column names."
            ),
            "Do not include hTERT, RIS, or CD columns in the RS comparison.",
            (
                "Use RS_TC_raw_TPM for the main TPM calculation. "
                "The project notes reserve RS_TC_batch_corrected_TPM "
                "for visualization."
            ),
            (
                "The R2 script must not run until the gene column and "
                "both sample groups have been completed and validated."
            ),
        ],
    }

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        json.dumps(config, indent=2) + "\n",
        encoding="utf-8",
    )

    logger.info("Created configuration template: %s", config_path)


def main() -> int:
    """Inspect the workbook, report columns, and prepare configuration."""
    args = parse_args()

    root = repo_root()
    input_path = resolve_from_repo(args.input, root)
    out_dir = resolve_from_repo(args.out_dir, root)

    inspection_dir = out_dir / "inspection"
    log_dir = out_dir / "logs"
    done_path = out_dir / "R1_rna_inspect.done.txt"
    config_path = out_dir / "rna_config.json"

    begin_run(done_path)

    logger, log_path = setup_logger(
        "R1_rna_inspect_tpm_table",
        log_dir,
    )

    if not input_path.is_file():
        logger.error("Workbook not found: %s", input_path)
        logger.error(
            "Place the workbook under 1_data_raw/ or supply --input."
        )
        return 2

    inspection_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Workbook: %s", input_path)
    logger.info("Inspection output: %s", inspection_dir)

    # Inspect workbook metadata and headers using Calamine.
    try:
        sheet_rows, column_rows, all_headers = read_workbook_headers(
            input_path,
            logger,
        )
    except Exception:
        logger.exception(
            "Could not inspect the workbook using the Calamine engine."
        )
        logger.error(
            "Check the workbook file and confirm that python-calamine "
            "is installed in this Python environment."
        )
        return 2

    # Save reports for every worksheet, even if the requested sheet is absent.
    pd.DataFrame(
        sheet_rows,
        columns=[
            "sheet_name",
            "reported_max_row",
            "reported_max_column",
            "header_count",
        ],
    ).to_csv(
        inspection_dir / "workbook_sheets.csv",
        index=False,
    )

    pd.DataFrame(
        column_rows,
        columns=[
            "sheet_name",
            "column_position_1based",
            "column_name",
        ],
    ).to_csv(
        inspection_dir / "columns_by_sheet.csv",
        index=False,
    )

    logger.info(
        "Saved worksheet report: %s",
        inspection_dir / "workbook_sheets.csv",
    )
    logger.info(
        "Saved column report: %s",
        inspection_dir / "columns_by_sheet.csv",
    )

    # Check that the requested worksheet exists.
    if args.sheet not in all_headers:
        logger.error("Sheet '%s' was not found.", args.sheet)
        logger.error(
            "Available sheets: %s",
            list(all_headers.keys()),
        )
        logger.error(
            "Inspect workbook_sheets.csv, then rerun with "
            "--sheet 'EXACT_NAME'."
        )
        return 2

    headers = all_headers[args.sheet]

    logger.info("Selected sheet: %s", args.sheet)
    logger.info("Header/column names (%d):", len(headers))

    for index, header in enumerate(headers, start=1):
        logger.info("  %03d | %s", index, header)

    # Identify potential sample columns without assigning sample groups.
    match_rows: list[dict] = []

    for header in headers:
        if not header:
            continue

        matched = [
            keyword
            for keyword in KEYWORDS
            if keyword.casefold() in header.casefold()
        ]

        if matched:
            match_rows.append(
                {
                    "column_name": header,
                    "matched_keywords": ";".join(matched),
                }
            )

    matches_path = inspection_dir / "sample_column_matches.csv"

    pd.DataFrame(
        match_rows,
        columns=["column_name", "matched_keywords"],
    ).to_csv(matches_path, index=False)

    logger.info("Keyword-matching columns (%d):", len(match_rows))

    for row in match_rows:
        logger.info(
            "  %s [matches: %s]",
            row["column_name"],
            row["matched_keywords"],
        )

    logger.info("Saved sample-column matches: %s", matches_path)

    # Create the config template only if it does not already exist.
    write_config_template(
        config_path=config_path,
        input_path=input_path,
        root=root,
        sheet_name=args.sheet,
        logger=logger,
    )

    write_done_marker(
        done_path,
        (
            "SUCCESS\n"
            f"WORKBOOK={input_path}\n"
            f"SHEET={args.sheet}\n"
            f"CONFIG={config_path}\n"
            f"LOG={log_path}\n"
        ),
    )

    logger.info("Wrote completion marker: %s", done_path)
    logger.info(
        "Next: check the printed headers and complete rna_config.json "
        "using confirmed sample columns."
    )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise

