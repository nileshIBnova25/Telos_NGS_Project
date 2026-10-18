from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "1_data_raw/5_motif_analysis"

df = pd.read_csv(RESULTS / "tead1_motif_region_scores.csv")
fasta = RESULTS / "benchmark_10000_sequences.fa"

# Read the extracted sequences.
sequences = {}
name = None
parts = []

with fasta.open() as handle:
    for line in handle:
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(parts).upper()
            name = line[1:].split("::", 1)[0].split("(", 1)[0]
            parts = []
        else:
            parts.append(line)

if name is not None:
    sequences[name] = "".join(parts).upper()

gc_values = []
valid_fractions = []

for _, row in df.iterrows():
    seq = sequences[row["region_id"]]
    valid = sum(base in "ACGT" for base in seq)
    gc = seq.count("G") + seq.count("C")

    gc_values.append(gc / valid if valid else float("nan"))
    valid_fractions.append(valid / len(seq) if seq else 0)

df["region_length"] = df["end"] - df["start"]
df["gc_fraction"] = gc_values
df["valid_fraction"] = valid_fractions

# Exclude only sequences with no valid DNA from motif comparisons.
usable = df[df["valid_fraction"] > 0].copy()

print("===== Sequence evaluability =====")
print(df.groupby("label")["valid_fraction"].agg(
    ["count", "mean", "min"]
).round(4))

print("\n===== Length and GC content by class =====")
print(usable.groupby("label").agg(
    n_regions=("region_id", "size"),
    mean_length=("region_length", "mean"),
    median_length=("region_length", "median"),
    mean_gc=("gc_fraction", "mean"),
    median_gc=("gc_fraction", "median"),
    motif_match_fraction=("has_motif_match", "mean"),
).round(4))

out = RESULTS / "tead1_motif_diagnostic_regions.csv"
df.to_csv(out, index=False)
print(f"\nSaved diagnostic table: {out}")
