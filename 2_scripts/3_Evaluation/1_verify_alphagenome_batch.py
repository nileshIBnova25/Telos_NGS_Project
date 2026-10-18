#!/usr/bin/env python3

"""Read-only validation of TELOS AlphaGenome batch outputs."""

import csv
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]

RUN_DIR = (
    PROJECT / "1_data_raw" / "3_alphagenome" / "02_atac_run"
)
BENCHMARK = (
    PROJECT / "1_data_raw" / "2_benchmark" / "benchmark_10000.bed"
)
PREDICTIONS = RUN_DIR / "benchmark_10000_alphagenome_atac_pred.csv"
MARKER = RUN_DIR / "step4B_alphagenome_atac_batch.done.txt"
LOG_DIR = RUN_DIR / "logs"


def count_bed_intervals(path):
    """Count non-empty, non-comment BED records."""
    if not path.exists():
        print(f"MISSING: {path}")
        return None

    with path.open() as handle:
        return sum(
            1 for line in handle
            if line.strip() and not line.startswith("#")
        )


def inspect_predictions(path):
    """Inspect CSV rows, columns, and blank values without pandas."""
    if not path.exists():
        print(f"MISSING: {path}")
        return None

    row_count = 0
    blank_rows = 0
    blank_values = {}
    columns = None

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames

        if not columns:
            print("ERROR: CSV header is missing.")
            return None

        blank_values = {column: 0 for column in columns}

        for row in reader:
            if not row or all(
                value is None or value.strip() == ""
                for value in row.values()
            ):
                blank_rows += 1
                continue

            row_count += 1

            for column in columns:
                value = row.get(column)
                if value is None or value.strip() == "":
                    blank_values[column] += 1

    print("CSV columns:", columns)
    print("CSV data rows:", row_count)
    print("Blank rows:", blank_rows)
    print("Blank values by column:", blank_values)

    return row_count


def main():
    print("TELOS AlphaGenome Step 4B — read-only audit")
    print("Project:", PROJECT)

    print("\n--- Benchmark ---")
    benchmark_count = count_bed_intervals(BENCHMARK)
    print("Benchmark intervals:", benchmark_count)

    print("\n--- Prediction CSV ---")
    prediction_count = inspect_predictions(PREDICTIONS)

    print("\n--- Completion marker ---")
    if MARKER.exists():
        print(MARKER.read_text().strip())
    else:
        print("Completion marker is missing.")

    print("\n--- Logs ---")
    if LOG_DIR.exists():
        logs = sorted(LOG_DIR.glob("step4B*.log"))
        if not logs:
            print("No Step 4B logs found.")
        for log in logs:
            print(f"\n[{log.name}]")
            lines = log.read_text(errors="replace").splitlines()
            for line in lines[-15:]:
                print(line)
    else:
        print("Log directory is missing.")

    print("\n--- Initial assessment ---")
    if benchmark_count is None or prediction_count is None:
        print("INCOMPLETE: Required benchmark or prediction file is missing.")
    elif prediction_count != benchmark_count:
        print("REVIEW REQUIRED: Prediction row count differs from benchmark.")
    else:
        print(
            "ROW COUNT MATCH ONLY: coordinate completeness, duplicates, "
            "prediction validity, and successful process completion "
            "still require verification."
        )


if __name__ == "__main__":
    main()
