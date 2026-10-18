from pathlib import Path
import sys

# Reuse the established calibration implementation, without modifying it.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import importlib.util
import pandas as pd

root = Path(__file__).resolve().parents[2]
source_script = (
    Path(__file__).resolve().parent / "17_calibrate_paper_guided_motifs.py"
)

# Load the established scanner/shuffle functions without running its main workflow.
source = source_script.read_text()
marker = "sequences = read_fasta(FASTA)"
if marker not in source:
    raise RuntimeError("Could not identify reusable functions in script 17.")

function_source = source.split(marker)[0]
namespace = {"__file__": str(source_script)}
exec(compile(function_source, str(source_script), "exec"), namespace)

read_fasta = namespace["read_fasta"]
load_pwm = namespace["load_pwm"]
score_max = namespace["score_max"]
dinucleotide_shuffle = namespace["dinucleotide_shuffle"]

import random
import numpy as np

folder = root / "1_data_raw/5_motif_analysis"
fasta = folder / "benchmark_10000_sequences.fa"
bed = root / "1_data_raw/2_benchmark/benchmark_10000.bed"
motif_dir = folder / "motifs/paper_guided_panel"
outdir = root / "3_results/tables/paper_guided_calibration"
outdir.mkdir(parents=True, exist_ok=True)

motifs = ["E2F1", "JUN", "NFKB1", "SMAD3"]
n_shuffles = 1000
rng = random.Random(20261017)

sequences = read_fasta(fasta)
labels = {}

with bed.open() as handle:
    for line in handle:
        if line.strip() and not line.startswith("#"):
            fields = line.rstrip("\n").split("\t")
            labels[fields[3]] = fields[4]

if set(sequences) != set(labels):
    raise ValueError("FASTA and BED region IDs do not match.")

all_rows = []

for motif in motifs:
    matrix_files = list(motif_dir.glob(f"{motif}_*.jaspar"))
    if len(matrix_files) != 1:
        raise ValueError(f"Expected one matrix for {motif}: {matrix_files}")

    pwm, motif_len = load_pwm(matrix_files[0])
    print(
        f"\nStarting {motif}: length={motif_len}, "
        f"shuffles={n_shuffles}",
        flush=True,
    )

    for i, (region_id, seq) in enumerate(sequences.items(), 1):
        if not seq or not set(seq) <= set("ACGT"):
            continue

        observed = score_max(seq, pwm, motif_len)
        null_scores = []

        for _ in range(n_shuffles):
            shuffled = dinucleotide_shuffle(seq, rng)
            if shuffled is not None:
                score = score_max(shuffled, pwm, motif_len)
                if score is not None:
                    null_scores.append(score)

        if observed is None or not null_scores:
            continue

        exceed = sum(value >= observed for value in null_scores)
        p = (exceed + 1) / (len(null_scores) + 1)

        all_rows.append({
            "motif": motif,
            "region_id": region_id,
            "label": labels[region_id],
            "sequence_length": len(seq),
            "observed_max_score": observed,
            "n_shuffles": len(null_scores),
            "null_mean_max_score": float(np.mean(null_scores)),
            "null_95pct_max_score": float(np.quantile(null_scores, 0.95)),
            "empirical_p": p,
        })

        if i % 500 == 0:
            print(f"{motif}: processed {i}/10000", flush=True)

    print(f"Finished {motif}", flush=True)

results = pd.DataFrame(all_rows)

def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    n = len(values)
    order = np.argsort(values)
    sorted_p = values[order]
    adjusted = np.minimum.accumulate(
        (sorted_p * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    output = np.empty(n)
    output[order] = adjusted
    return output

results["BH_q_value_within_motif"] = np.nan
for motif in motifs:
    mask = results["motif"] == motif
    results.loc[mask, "BH_q_value_within_motif"] = bh_adjust(
        results.loc[mask, "empirical_p"].to_numpy()
    )

outfile = outdir / "remaining_four_motif_region_calibration.csv"
results.to_csv(outfile, index=False)

summary = results.groupby(["motif", "label"]).agg(
    n_regions=("region_id", "size"),
    median_observed_score=("observed_max_score", "median"),
    median_null_mean_score=("null_mean_max_score", "median"),
    n_empirical_p_le_005=("empirical_p", lambda x: int((x <= 0.05).sum())),
    n_BH_q_le_005=("BH_q_value_within_motif", lambda x: int((x <= 0.05).sum())),
).reset_index()

summary_file = outdir / "remaining_four_motif_class_summary.csv"
summary.to_csv(summary_file, index=False)

print("\nClass summary:")
print(summary.round(4).to_string(index=False))
print("\nSaved:", outfile)
print("Saved:", summary_file)
