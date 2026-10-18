#!/usr/bin/env python3

import csv
import math
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BED = ROOT / "1_data_raw/2_benchmark/benchmark_10000.bed"
FASTA = ROOT / "1_data_raw/0_reference/hg38/hg38.fa"
MOTIF = ROOT / "1_data_raw/5_motif_analysis/motifs/TEAD1_MA0090.4.jaspar"
OUT = ROOT / "1_data_raw/5_motif_analysis"
OUT.mkdir(parents=True, exist_ok=True)

# Load the JASPAR position-frequency matrix.
matrix = {}
with MOTIF.open() as handle:
    for line in handle:
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        base, values = line.split("[", 1)
        matrix[base.strip()] = np.array(
            [float(x) for x in values.rstrip(" ]").split()],
            dtype=float,
        )

if set(matrix) != {"A", "C", "G", "T"}:
    raise ValueError("Expected A, C, G and T rows in the JASPAR matrix.")

lengths = {len(values) for values in matrix.values()}
if len(lengths) != 1:
    raise ValueError("The JASPAR matrix rows have inconsistent lengths.")

motif_length = lengths.pop()
if motif_length == 0:
    raise ValueError("The motif matrix is empty.")

# Convert counts to log2 odds, normalising across bases at each position.
pseudocount = 0.5
bases = "ACGT"
count_matrix = np.vstack([matrix[base] for base in bases])
position_totals = count_matrix.sum(axis=0) + 4 * pseudocount

probabilities = {
    base: (matrix[base] + pseudocount) / position_totals
    for base in bases
}

log_odds = {
    base: np.log2(probabilities[base] / 0.25)
    for base in bases
}

# The reverse-complement score is evaluated explicitly.
complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

def score_sequence(sequence):
    sequence = sequence.upper()
    if len(sequence) != motif_length or any(
        base not in "ACGT" for base in sequence
    ):
        return None

    forward = sum(
        log_odds[base][i] for i, base in enumerate(sequence)
    )
    reverse_sequence = "".join(complement[b] for b in sequence[::-1])
    reverse = sum(
        log_odds[base][i] for i, base in enumerate(reverse_sequence)
    )
    return max(float(forward), float(reverse))

def read_bed(path):
    records = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"Expected at least 5 BED columns: {line}")
            chrom, start, end, region_id, label = fields[:5]
            records.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "region_id": region_id,
                "label": label,
            })
    return records

records = read_bed(BED)
if len(records) != 10000:
    raise ValueError(f"Expected 10,000 BED intervals; found {len(records)}")

# Extract the benchmark sequences from the local reference with bedtools.
# Use a temporary BED in the output directory and preserve BED coordinates.
sequence_bed = OUT / "benchmark_10000_for_fasta.bed"
sequence_fasta = OUT / "benchmark_10000_sequences.fa"

with sequence_bed.open("w") as handle:
    for rec in records:
        handle.write(
            f"{rec['chrom']}\t{rec['start']}\t{rec['end']}\t"
            f"{rec['region_id']}\n"
        )

try:
    with sequence_fasta.open("w") as handle:
        subprocess.run(
            ["bedtools", "getfasta", "-fi", str(FASTA),
             "-bed", str(sequence_bed), "-name"],
            check=True,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
        )
except subprocess.CalledProcessError as exc:
    raise RuntimeError(
        "bedtools getfasta failed. Check chromosome naming and hg38 reference."
    ) from exc

sequences = {}
with sequence_fasta.open() as handle:
    current_id = None
    chunks = []
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if current_id is not None:
                sequences[current_id] = "".join(chunks).upper()
            current_id = line[1:].split("::", 1)[0].split("(", 1)[0]
            chunks = []
        else:
            chunks.append(line)
    if current_id is not None:
        sequences[current_id] = "".join(chunks).upper()

# Scan every possible motif-length window. A site is counted if its score
# reaches the chosen threshold; we use a transparent exploratory threshold
# of 80% of the maximum possible log-odds score.
max_score = sum(
    max(log_odds[base][i] for base in "ACGT")
    for i in range(motif_length)
)
threshold = 0.8 * max_score

results = []
for rec in records:
    sequence = sequences.get(rec["region_id"])
    if sequence is None:
        raise ValueError(f"Missing extracted sequence: {rec['region_id']}")

    scores = []
    for i in range(max(0, len(sequence) - motif_length + 1)):
        score = score_sequence(sequence[i:i + motif_length])
        if score is not None:
            scores.append(score)

    matches = [score for score in scores if score >= threshold]
    results.append({
        **rec,
        "sequence_length": len(sequence),
        "valid_windows": len(scores),
        "max_motif_score": max(scores) if scores else "",
        "motif_match_count": len(matches),
        "has_motif_match": int(bool(matches)),
    })

out_csv = OUT / "tead1_motif_region_scores.csv"
with out_csv.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(results[0].keys()))
    writer.writeheader()
    writer.writerows(results)

summary_csv = OUT / "tead1_motif_class_summary.csv"
labels = sorted({rec["label"] for rec in results})
with summary_csv.open("w", newline="") as handle:
    fields = [
        "label", "n_regions", "regions_with_match",
        "fraction_with_match", "mean_match_count", "median_max_score",
    ]
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    for label in labels:
        group = [r for r in results if r["label"] == label]
        writer.writerow({
            "label": label,
            "n_regions": len(group),
            "regions_with_match": sum(r["has_motif_match"] for r in group),
            "fraction_with_match": np.mean(
                [r["has_motif_match"] for r in group]
            ),
            "mean_match_count": np.mean(
                [r["motif_match_count"] for r in group]
            ),
            "median_max_score": np.median([
                r["max_motif_score"] for r in group
                if r["max_motif_score"] != ""
            ]),
        })

print(f"Motif length: {motif_length}")
print(f"Exploratory score threshold: {threshold:.4f} bits")
print(f"Regions scanned: {len(results)}")
print(f"Region-level output: {out_csv}")
print(f"Class summary: {summary_csv}")
print("\nClass summary:")
for label in labels:
    group = [r for r in results if r["label"] == label]
    n_match = sum(r["has_motif_match"] for r in group)
    print(f"{label}: {n_match}/{len(group)} regions with a motif match")
