from pathlib import Path
from collections import Counter
import csv
import random

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "1_data_raw/5_motif_analysis"
FASTA = FOLDER / "benchmark_10000_sequences.fa"
MOTIF = FOLDER / "motifs/TEAD1_MA0090.4.jaspar"

N_SHUFFLES = 1000
SEED = 20261017
PSEUDOCOUNT = 0.5
BATCH_SIZE = 100

# Load and validate the JASPAR matrix.
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

bases = "ACGT"
if set(matrix) != set(bases):
    raise ValueError("Expected A, C, G and T in the PWM.")

count_matrix = np.vstack([matrix[b] for b in bases])
totals = count_matrix.sum(axis=0) + 4 * PSEUDOCOUNT
prob = {
    b: (matrix[b] + PSEUDOCOUNT) / totals
    for b in bases
}
log_odds = {
    b: np.log2(prob[b] / 0.25)
    for b in bases
}
MOTIF_LEN = len(matrix["A"])
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}

def score_window(seq):
    if len(seq) != MOTIF_LEN or any(b not in "ACGT" for b in seq):
        return None
    forward = sum(log_odds[b][i] for i, b in enumerate(seq))
    reverse_seq = "".join(COMPLEMENT[b] for b in seq[::-1])
    reverse = sum(log_odds[b][i] for i, b in enumerate(reverse_seq))
    return max(float(forward), float(reverse))

def max_score(seq):
    best = None
    for i in range(len(seq) - MOTIF_LEN + 1):
        score = score_window(seq[i:i + MOTIF_LEN])
        if score is not None and (best is None or score > best):
            best = score
    return best

def dinucleotide_counts(seq):
    return Counter(seq[i:i+2] for i in range(len(seq)-1))

def dinucleotide_shuffle(seq, rng):
    # Randomized Eulerian trail: preserves all dinucleotide counts.
    edges = {}
    for a, b in zip(seq[:-1], seq[1:]):
        edges.setdefault(a, []).append(b)
    for a in edges:
        rng.shuffle(edges[a])

    stack = [seq[0]]
    trail = []
    while stack:
        a = stack[-1]
        if edges.get(a):
            stack.append(edges[a].pop())
        else:
            trail.append(stack.pop())

    shuffled = "".join(reversed(trail))
    if len(shuffled) != len(seq):
        raise RuntimeError("Shuffled sequence length changed.")
    if dinucleotide_counts(shuffled) != dinucleotide_counts(seq):
        raise RuntimeError("Dinucleotide counts were not preserved.")
    return shuffled

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

sequences = read_fasta(FASTA)
scores_file = FOLDER / "tead1_motif_region_scores.csv"
original = pd.read_csv(scores_file)

if len(sequences) != 10000 or len(original) != 10000:
    raise ValueError("Expected exactly 10,000 sequences and score rows.")

if original["region_id"].duplicated().any():
    raise ValueError("Duplicate region IDs in original score table.")

if set(sequences) != set(original["region_id"]):
    raise ValueError("FASTA IDs and score-table IDs do not match.")

rng = random.Random(SEED)
rows = []

for start in range(0, len(original), BATCH_SIZE):
    batch = original.iloc[start:start+BATCH_SIZE]
    for _, rec in batch.iterrows():
        region_id = rec["region_id"]
        seq = sequences[region_id]

        # Regions containing N or other ambiguous bases are unscorable.
        if not seq or any(b not in "ACGT" for b in seq):
            rows.append({
                "region_id": region_id,
                "label": rec["label"],
                "sequence_length": len(seq),
                "observed_max_score": np.nan,
                "n_shuffles": 0,
                "null_mean_max_score": np.nan,
                "null_95pct_max_score": np.nan,
                "empirical_p": np.nan,
            })
            continue

        observed = max_score(seq)
        null_scores = []

        for _ in range(N_SHUFFLES):
            shuffled = dinucleotide_shuffle(seq, rng)
            score = max_score(shuffled)
            if score is not None:
                null_scores.append(score)

        if observed is None or not null_scores:
            p_empirical = np.nan
            null_mean = np.nan
            null_95 = np.nan
        else:
            # One-sided test: unusually high observed maximum score.
            exceed = sum(s >= observed for s in null_scores)
            p_empirical = (exceed + 1) / (len(null_scores) + 1)
            null_mean = float(np.mean(null_scores))
            null_95 = float(np.quantile(null_scores, 0.95))

        rows.append({
            "region_id": region_id,
            "label": rec["label"],
            "sequence_length": len(seq),
            "observed_max_score": observed,
            "n_shuffles": len(null_scores),
            "null_mean_max_score": null_mean,
            "null_95pct_max_score": null_95,
            "empirical_p": p_empirical,
        })

    done = min(start + BATCH_SIZE, len(original))
    print(f"Processed {done}/{len(original)} regions", flush=True)

out = pd.DataFrame(rows)

# Benjamini-Hochberg adjustment over scorable regions.
valid = out["empirical_p"].notna()
pvals = out.loc[valid, "empirical_p"].to_numpy()
order = np.argsort(pvals)
ranked = pvals[order]
m = len(ranked)

adjusted = np.minimum.accumulate(
    (ranked * m / np.arange(1, m + 1))[::-1]
)[::-1]
adjusted = np.clip(adjusted, 0, 1)

qvals = np.empty_like(adjusted)
qvals[order] = adjusted
out["BH_q_value"] = np.nan
out.loc[valid, "BH_q_value"] = qvals

OUTDIR = ROOT / "3_results/tables/tead1_calibration_1000"
OUTDIR.mkdir(parents=True, exist_ok=True)
region_output = OUTDIR / "tead1_dinucleotide_calibration_regions_1000.csv"
out.to_csv(region_output, index=False)

summary = (
    out.groupby("label", dropna=False)
    .agg(
        n_regions=("region_id", "size"),
        n_scorable=("empirical_p", "count"),
        mean_observed_score=("observed_max_score", "mean"),
        mean_null_score=("null_mean_max_score", "mean"),
        fraction_p_le_005=("empirical_p", lambda x: (x <= 0.05).mean()),
        fraction_q_le_005=("BH_q_value", lambda x: (x <= 0.05).mean()),
    )
    .reset_index()
)

summary_output = OUTDIR / "tead1_dinucleotide_calibration_summary_1000.csv"
summary.to_csv(summary_output, index=False)

print("\nCalibration summary:")
print(summary.round(4).to_string(index=False))
print("\nSaved region results:", region_output)
print("Saved class summary:", summary_output)
print(
    "\nInterpretation note: empirical p-values have resolution limited by "
    f"{N_SHUFFLES} shuffles. Treat this run as exploratory, not definitive."
)
