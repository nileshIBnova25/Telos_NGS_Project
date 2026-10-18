from pathlib import Path
from collections import Counter
import random

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "1_data_raw/5_motif_analysis"
FASTA = FOLDER / "benchmark_10000_sequences.fa"
MOTIF_DIR = FOLDER / "motifs/paper_guided_panel"
OUTDIR = ROOT / "3_results/tables/paper_guided_calibration"
OUTDIR.mkdir(parents=True, exist_ok=True)

MOTIFS = ["CEBPA", "FOS", "RELA", "SMAD4"]
N_SHUFFLES = 1000
SEED = 20261016
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

def dinucleotide_counts(seq):
    return Counter(seq[i:i+2] for i in range(len(seq)-1))

def dinucleotide_shuffle(seq, rng):
    if len(seq) < 2 or not set(seq) <= set(BASES):
        return None

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
        raise ValueError(f"Invalid matrix: {path}")

    lengths = {len(counts[b]) for b in BASES}
    if len(lengths) != 1:
        raise ValueError(f"Matrix rows have different lengths: {path}")

    matrix = np.vstack([counts[b] for b in BASES])
    totals = matrix.sum(axis=0) + 4 * PSEUDOCOUNT
    probabilities = {
        b: (counts[b] + PSEUDOCOUNT) / totals for b in BASES
    }
    log_odds = {
        b: np.log2(probabilities[b] / 0.25) for b in BASES
    }
    motif_len = lengths.pop()
    return log_odds, motif_len

def score_max(seq, log_odds, motif_len):
    best = None
    for i in range(len(seq) - motif_len + 1):
        window = seq[i:i + motif_len]
        if any(b not in BASES for b in window):
            continue

        forward = sum(log_odds[b][j] for j, b in enumerate(window))
        reverse_seq = "".join(COMPLEMENT[b] for b in window[::-1])
        reverse = sum(
            log_odds[b][j] for j, b in enumerate(reverse_seq)
        )
        score = max(float(forward), float(reverse))

        if best is None or score > best:
            best = score
    return best

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
bed = ROOT / "1_data_raw/2_benchmark/benchmark_10000.bed"
metadata = {}

with bed.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        metadata[fields[3]] = fields[4]

if len(sequences) != 10000 or len(metadata) != 10000:
    raise ValueError("Expected 10,000 sequences and 10,000 BED records.")
if set(sequences) != set(metadata):
    raise ValueError("FASTA and BED region IDs do not match.")

rng = random.Random(SEED)
all_results = []

for motif in MOTIFS:
    matrix_files = list(MOTIF_DIR.glob(f"{motif}_*.jaspar"))
    if len(matrix_files) != 1:
        raise ValueError(f"Expected one matrix for {motif}; found {matrix_files}")

    log_odds, motif_len = load_pwm(matrix_files[0])
    print(
        f"\nStarting {motif}; motif length={motif_len}; "
        f"shuffles per region={N_SHUFFLES}",
        flush=True,
    )

    for index, (region_id, seq) in enumerate(sequences.items(), start=1):
        label = metadata[region_id]

        if not seq or not set(seq) <= set(BASES):
            all_results.append({
                "motif": motif,
                "region_id": region_id,
                "label": label,
                "sequence_length": len(seq),
                "observed_max_score": np.nan,
                "n_shuffles": 0,
                "null_mean_max_score": np.nan,
                "null_95pct_max_score": np.nan,
                "empirical_p": np.nan,
            })
            continue

        observed = score_max(seq, log_odds, motif_len)
        null_scores = []

        for _ in range(N_SHUFFLES):
            shuffled = dinucleotide_shuffle(seq, rng)
            if shuffled is not None:
                score = score_max(shuffled, log_odds, motif_len)
                if score is not None:
                    null_scores.append(score)

        if observed is None or not null_scores:
            pvalue = np.nan
            null_mean = np.nan
            null_95 = np.nan
        else:
            exceed = sum(s >= observed for s in null_scores)
            pvalue = (exceed + 1) / (len(null_scores) + 1)
            null_mean = float(np.mean(null_scores))
            null_95 = float(np.quantile(null_scores, 0.95))

        all_results.append({
            "motif": motif,
            "region_id": region_id,
            "label": label,
            "sequence_length": len(seq),
            "observed_max_score": observed,
            "n_shuffles": len(null_scores),
            "null_mean_max_score": null_mean,
            "null_95pct_max_score": null_95,
            "empirical_p": pvalue,
        })

        if index % 500 == 0:
            print(f"{motif}: processed {index}/10000", flush=True)

    print(f"Finished {motif}", flush=True)

results = pd.DataFrame(all_results)

# BH correction separately for each motif across its scorable regions.
results["BH_q_value_within_motif"] = np.nan
for motif in MOTIFS:
    mask = (
        (results["motif"] == motif)
        & results["empirical_p"].notna()
    )
    if mask.any():
        results.loc[mask, "BH_q_value_within_motif"] = bh_adjust(
            results.loc[mask, "empirical_p"].to_numpy()
        )

region_path = OUTDIR / "four_motif_region_calibration.csv"
results.to_csv(region_path, index=False)

summary_rows = []
for motif in MOTIFS:
    subset = results[results["motif"] == motif]
    for label, group in subset.groupby("label"):
        valid = group["empirical_p"].notna()
        summary_rows.append({
            "motif": motif,
            "label": label,
            "n_regions": len(group),
            "n_scorable": int(valid.sum()),
            "median_observed_max_score": group["observed_max_score"].median(),
            "median_null_mean_max_score": group["null_mean_max_score"].median(),
            "n_empirical_p_le_0.05": int((group.loc[valid, "empirical_p"] <= 0.05).sum()),
            "n_BH_q_le_0.05": int((group["BH_q_value_within_motif"] <= 0.05).sum()),
        })

summary = pd.DataFrame(summary_rows)
summary_path = OUTDIR / "four_motif_class_calibration_summary.csv"
summary.to_csv(summary_path, index=False)

print("\nClass calibration summary:")
print(summary.round(4).to_string(index=False))
print("\nSaved region results:", region_path)
print("Saved class summary:", summary_path)
print(
    "\nImportant: empirical p-values have minimum resolution "
    f"1/{N_SHUFFLES + 1}. BH correction is within each motif across "
    "regions. This calibration tests unusual maximum scores conditional "
    "on each sequence's dinucleotide composition; it does not prove "
    "transcription-factor binding."
)
