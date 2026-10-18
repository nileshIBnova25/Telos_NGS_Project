from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "1_data_raw/5_motif_analysis"
FASTA = FOLDER / "benchmark_10000_sequences.fa"
BED = ROOT / "1_data_raw/2_benchmark/benchmark_10000.bed"
SCORES = ROOT / "3_results/tables/paper_guided_motif_region_scores.csv"
OUT = ROOT / "3_results/tables/paper_guided_motif_adjusted_tests.csv"

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
metadata = []

with BED.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        region_id, label = fields[3], fields[4]
        seq = sequences[region_id]
        valid_bases = sum(b in "ACGT" for b in seq)
        gc_bases = sum(b in "GC" for b in seq)

        metadata.append({
            "region_id": region_id,
            "label": label,
            "length": len(seq),
            "log_length": np.log(len(seq)) if seq else np.nan,
            "GC_percent": 100 * gc_bases / valid_bases
                if valid_bases else np.nan,
            "valid_base_fraction": valid_bases / len(seq)
                if seq else np.nan,
        })

meta = pd.DataFrame(metadata)
scores = pd.read_csv(SCORES)
df = scores.merge(meta, on=["region_id", "label"], validate="many_to_one")

if len(df) != 80000:
    raise ValueError(f"Expected 80,000 motif-region rows; found {len(df)}")

# Exclude regions with ambiguous bases from this adjusted comparison.
df = df[df["valid_base_fraction"] == 1.0].copy()
rows = []

for motif in sorted(df["motif"].unique()):
    for target in ["RS_gained", "RS_lost"]:
        subset = df[df["label"].isin([target, "Stable"])].copy()
        subset["is_target"] = (subset["label"] == target).astype(int)
        subset["log_length_z"] = (
            subset["log_length"] - subset["log_length"].mean()
        ) / subset["log_length"].std()
        subset["GC_z"] = (
            subset["GC_percent"] - subset["GC_percent"].mean()
        ) / subset["GC_percent"].std()

        motif_data = subset[subset["motif"] == motif].copy()
        X = sm.add_constant(
            motif_data[["is_target", "log_length_z", "GC_z"]]
        )
        y = motif_data["has_match"]

        try:
            fit = sm.GLM(
                y, X, family=sm.families.Binomial()
            ).fit()

            beta = fit.params["is_target"]
            se = fit.bse["is_target"]
            pvalue = fit.pvalues["is_target"]

            rows.append({
                "motif": motif,
                "comparison": f"{target}_vs_Stable",
                "n_regions": len(motif_data),
                "target_match_fraction": motif_data.loc[
                    motif_data["is_target"] == 1, "has_match"
                ].mean(),
                "stable_match_fraction": motif_data.loc[
                    motif_data["is_target"] == 0, "has_match"
                ].mean(),
                "adjusted_odds_ratio": np.exp(beta),
                "CI95_lower": np.exp(beta - 1.96 * se),
                "CI95_upper": np.exp(beta + 1.96 * se),
                "p_value": pvalue,
                "status": "OK",
            })
        except Exception as exc:
            rows.append({
                "motif": motif,
                "comparison": f"{target}_vs_Stable",
                "n_regions": len(motif_data),
                "status": f"FAILED: {exc}",
            })

result = pd.DataFrame(rows)
valid = result["status"] == "OK"
result["BH_q_value"] = np.nan

if valid.any():
    p = result.loc[valid, "p_value"].to_numpy()
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    q = np.empty(n)
    q[order] = adjusted
    result.loc[valid, "BH_q_value"] = q

result.to_csv(OUT, index=False)

print(result.round(5).to_string(index=False))
print("\nSaved:", OUT)
print(
    "\nInterpretation: odds ratios are adjusted for log sequence length "
    "and GC content, among fully unambiguous sequences. These tests still "
    "depend on the exploratory 80%-of-maximum motif threshold."
)
