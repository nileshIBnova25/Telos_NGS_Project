#!/usr/bin/env python3

"""Read-only coordinate validation for TELOS AlphaGenome predictions."""

import csv
from collections import Counter
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]

BENCHMARK = (
    PROJECT / "1_data_raw" / "2_benchmark" / "benchmark_10000.bed"
)

PREDICTIONS = (
    PROJECT / "1_data_raw" / "3_alphagenome" / "02_atac_run"
    / "benchmark_10000_alphagenome_atac_pred.csv"
)


def read_bed(path):
    """Read chromosome, start, and end from a BED file."""
    intervals = []

    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 3:
                raise ValueError(
                    f"Invalid BED record at line {line_number}: "
                    "fewer than 3 tab-separated columns."
                )

            try:
                chrom = fields[0].strip()
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"Invalid coordinates at BED line {line_number}."
                ) from error

            if not chrom or start < 0 or end <= start:
                raise ValueError(
                    f"Invalid interval at BED line {line_number}: {line}"
                )

            intervals.append((chrom, start, end))

    return intervals


def read_predictions(path):
    """Read and validate genomic intervals from the prediction CSV."""
    intervals = []
    malformed = 0
    missing_signal = 0
    signal_errors = 0
    row_count = 0

    required = {"chr", "start", "end", "pred_atac_mean"}

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)

        if not reader.fieldnames or not required.issubset(
            set(reader.fieldnames)
        ):
            raise ValueError(
                "Prediction CSV is missing required columns. "
                f"Found: {reader.fieldnames}"
            )

        for row_number, row in enumerate(reader, start=2):
            row_count += 1

            try:
                chrom = row["chr"].strip()
                start = int(row["start"])
                end = int(row["end"])

                if not chrom or start < 0 or end <= start:
                    raise ValueError("Invalid interval")

                intervals.append((chrom, start, end))

            except (ValueError, AttributeError, TypeError):
                malformed += 1
                continue

            signal = row.get("pred_atac_mean")

            if signal is None or not signal.strip():
                missing_signal += 1
            else:
                try:
                    value = float(signal)
                    if not __import__("math").isfinite(value):
                        raise ValueError("Non-finite signal")
                except ValueError:
                    signal_errors += 1

    return (
        intervals,
        row_count,
        malformed,
        missing_signal,
        signal_errors,
    )


def duplicate_count(intervals):
    counts = Counter(intervals)
    return sum(count - 1 for count in counts.values() if count > 1)


def main():
    print("TELOS — AlphaGenome coordinate validation")
    print("Mode: read-only; no source files will be modified.")

    for path in (BENCHMARK, PREDICTIONS):
        if not path.is_file():
            print(f"\nERROR: Required file is missing: {path}")
            return

    try:
        benchmark = read_bed(BENCHMARK)
        (
            predicted,
            csv_rows,
            malformed,
            missing_signal,
            signal_errors,
        ) = read_predictions(PREDICTIONS)
    except (OSError, ValueError) as error:
        print(f"\nVALIDATION ERROR: {error}")
        return

    benchmark_set = set(benchmark)
    predicted_set = set(predicted)

    missing = benchmark_set - predicted_set
    unexpected = predicted_set - benchmark_set

    print("\n--- Record counts ---")
    print("Benchmark BED records:", len(benchmark))
    print("Unique benchmark intervals:", len(benchmark_set))
    print("Prediction CSV data rows:", csv_rows)
    print("Valid prediction intervals:", len(predicted))
    print("Unique prediction intervals:", len(predicted_set))

    print("\n--- Coordinate checks ---")
    print("Duplicate benchmark records:", duplicate_count(benchmark))
    print("Duplicate prediction records:", duplicate_count(predicted))
    print("Missing benchmark intervals:", len(missing))
    print("Unexpected prediction intervals:", len(unexpected))
    print("Malformed prediction rows:", malformed)

    print("\n--- Prediction signal checks ---")
    print("Missing prediction signals:", missing_signal)
    print("Invalid or non-finite signals:", signal_errors)

    if missing:
        print("\nExamples of missing intervals:")
        for interval in sorted(missing)[:10]:
            print(interval)

    if unexpected:
        print("\nExamples of unexpected intervals:")
        for interval in sorted(unexpected)[:10]:
            print(interval)

    passed = (
        len(benchmark) == 10000
        and len(benchmark_set) == 10000
        and csv_rows == 10000
        and len(predicted) == 10000
        and len(predicted_set) == 10000
        and duplicate_count(predicted) == 0
        and duplicate_count(benchmark) == 0
        and not missing
        and not unexpected
        and malformed == 0
        and missing_signal == 0
        and signal_errors == 0
    )

    print("\n--- Overall assessment ---")
    if passed:
        print("PASS: All 10,000 benchmark intervals match uniquely.")
    else:
        print("REVIEW REQUIRED: One or more validation checks failed.")
        print("No source files were modified.")


if __name__ == "__main__":
    main()
