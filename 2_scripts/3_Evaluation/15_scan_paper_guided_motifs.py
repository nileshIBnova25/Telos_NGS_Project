#!/usr/bin/env python3

from pathlib import Path
import csv
import math

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "1_data_raw/5_motif_analysis"
MOTIF_DIR = FOLDER / "motifs/paper_guided_panel"
FASTA = FOLDER / "benchmark_10000_sequences.fa"
OUTDIR = ROOT / "3_results/tables"
OUTDIR.mkdir(parents=True, exist_ok=True)

PSEUDOCOUNT = 0.5
BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

def read_fasta(path):
    sequences = {}
    name = None
    chunks = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks).upper()
                name = line[1:].split("::", 1)[0].split("(", 1)[0]
                chunks = []
            else:
                chunks.append(line)

        if name is not None:
            sequences[name] = "".join(chunks).upper()

    return sequences

def read_bed(path):
    records = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            records.append({
                "region_id": fields[3],
                "label": fields[4],
            })
    return records

def load_pwm(path):
    counts = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            base, values = line.split("[", 1)
            counts[base.strip()] = np.array(
                [float(x) for x in values.rstrip(" ]").split()],
                dtype=float,
            )

    if set(counts) != set(BASES):
        raise ValueError(f"Invalid JASPAR matrix: {path}")

    lengths = {len(counts[b]) for b in BASES}
    if len(lengths) != 1:
        raise ValueError(f"Inconsistent matrix lengths: {path}")

    count_matrix = np.vstack([counts[b] for b in BASES])
    totals = count_matrix.sum(axis=0) + 4 * PSEUDOCOUNT
    probabilities = {
        b: (counts[b] + PSEUDOCOUNT) / totals
        for b in BASES
    }
    log_odds = {
        b: np.log2(probabilities[b] / 0.25)
        for b in BASES
    }

    motif_len = lengths.pop()
    max_score = sum(
        max(log_odds[b][i] for b in BASES)
        for i in range(motif_len)
    )

    return log_odds, motif_len, 0.8 * max_score

def scan_sequence(seq, log_odds, motif_len):
    best = None
    n_matches = 0
    threshold = 0.8 * sum(
        max(log_odds[b][i] for b in BASES)
        for i in range(motif_len)
    )

    for start in range(len(seq) - motif_len + 1):
        window = seq[start:start + motif_len]
        if any(b not in BASES for b in window):
            continue

        forward = sum(log_odds[b][i] for i, b in enumerate(window))
        reverse_seq = "".join(COMPLEMENT[b] for b in window[::-1])
        reverse = sum(
            log_odds[b][i] for i, b in enumerate(reverse_seq)
        )
        score = max(float(forward), float(reverse))

        if best is None or score > best:
            best = score
        if score >= threshold:
            n_matches += 1

    return best, n_matches, threshold

def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty(n)
    result[order] = adjusted
    return result

sequences = read_fasta(FASTA)
records = read_bed(ROOT / "1_data_raw/2_benchmark/benchmark_10000.bed")

if len(sequences) != 10000 or len(records) != 10000:
    raise ValueError(
        f"Expected 10,000 sequences and BED records; "
        f"got {len(sequences)} and {len(records)}"
    )

if set(sequences) != {r["region_id"] for r in records}:
    raise ValueError("Sequence IDs and BED region IDs do not match.")

motif_files = sorted(MOTIF_DIR.glob("*.jaspar"))
if len(motif_files) != 8:
    raise ValueError(
        f"Expected 8 motif files in {MOTIF_DIR}; found {len(motif_files)}"
    )

region_rows = []
threshold_rows = []

for motif_file in motif_files:
    motif_name = motif_file.name.split("_MA", 1)[0]
    log_odds, motif_len, threshold = load_pwm(motif_file)
    threshold_rows.append({
        "motif": motif_name,
        "matrix_file": motif_file.name,
        "motif_length": motif_len,
        "threshold_bits": threshold,
        "threshold_rule": "80% of maximum possible log-odds score",
    })

    for rec in records:
        seq = sequences[rec["region_id"]]
        best, n_matches, _ = scan_sequence(seq, log_odds, motif_len)
        region_rows.append({
            **rec,
            "motif": motif_name,
            "sequence_length": len(seq),
            "max_score_bits": best if best is not None else np.nan,
            "match_count": n_matches,
            "has_match": int(n_matches > 0),
        })

    print(
        f"Scanned {motif_name} ({motif_len} bp); "
        f"threshold={threshold:.3f} bits",
        flush=True,
    )

region_df = pd.DataFrame(region_rows)
region_path = OUTDIR / "paper_guided_motif_region_scores.csv"
region_df.to_csv(region_path, index=False)

pd.DataFrame(threshold_rows).to_csv(
    OUTDIR / "paper_guided_motif_thresholds.csv", index=False
)

class_summary = (
    region_df.groupby(["motif", "label"])
    .agg(
        n_regions=("region_id", "size"),
        regions_with_match=("has_match", "sum"),
        mean_match_count=("match_count", "mean"),
        median_max_score=("max_score_bits", "median"),
    )
    .reset_index()
)
class_summary["fraction_with_match"] = (
    class_summary["regions_with_match"] / class_summary["n_regions"]
)
class_summary.to_csv(
    OUTDIR / "paper_guided_motif_class_summary.csv", index=False
)

# Planned comparisons: gained vs Stable, and lost vs Stable.
tests = []
for motif in sorted(region_df["motif"].unique()):
    subset = region_df[region_df["motif"] == motif]

    for target in ["RS_gained", "RS_lost"]:
        a = subset[subset["label"] == target]["has_match"]
        b = subset[subset["label"] == "Stable"]["has_match"]

        if len(a) == 0 or len(b) == 0:
            continue

        a_yes, a_no = int(a.sum()), int(len(a) - a.sum())
        b_yes, b_no = int(b.sum()), int(len(b) - b.sum())

        odds_ratio, pvalue = fisher_exact(
            [[a_yes, a_no], [b_yes, b_no]],
            alternative="two-sided",
        )

        tests.append({
            "motif": motif,
            "comparison": f"{target}_vs_Stable",
            "target_positive": a_yes,
            "target_total": len(a),
            "target_fraction": a_yes / len(a),
            "stable_positive": b_yes,
            "stable_total": len(b),
            "stable_fraction": b_yes / len(b),
            "odds_ratio": odds_ratio,
            "p_value": pvalue,
        })

test_df = pd.DataFrame(tests)
test_df["BH_q_value"] = bh_adjust(test_df["p_value"].to_numpy())
test_df.to_csv(
    OUTDIR / "paper_guided_motif_class_tests.csv", index=False
)

print("\nClass frequencies:")
print(class_summary.round(4).to_string(index=False))
print("\nPlanned class comparisons:")
print(test_df.round(5).to_string(index=False))
print("\nSaved outputs:")
print(region_path)
print(OUTDIR / "paper_guided_motif_thresholds.csv")
print(OUTDIR / "paper_guided_motif_class_summary.csv")
print(OUTDIR / "paper_guided_motif_class_tests.csv")
print(
    "\nCaution: 80%-of-maximum thresholds are exploratory and do not "
    "equalize false-positive rates across motifs. Class-test q-values "
    "are conditional on this threshold choice."
)
