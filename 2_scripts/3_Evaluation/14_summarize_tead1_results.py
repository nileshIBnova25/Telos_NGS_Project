from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "1_data_raw/5_motif_analysis"

calibration = pd.read_csv(
    FOLDER / "tead1_dinucleotide_calibration_summary.csv"
)
tests = pd.read_csv(
    FOLDER / "tead1_class_difference_tests.csv"
)
diagnostics = pd.read_csv(
    FOLDER / "tead1_motif_diagnostic_regions.csv"
)

diagnostics = diagnostics.loc[
    diagnostics["valid_fraction"] > 0
].copy()

diagnostics["has_motif_match"] = (
    diagnostics["has_motif_match"].astype(int)
)

match_summary = (
    diagnostics.groupby("label")["has_motif_match"]
    .agg(n_evaluable="count", n_with_match="sum", match_fraction="mean")
    .reset_index()
)

summary = calibration.merge(match_summary, on="label", how="left")

output = FOLDER / "tead1_final_class_summary.csv"
summary.to_csv(output, index=False)

print("TEAD1 motif analysis: consolidated class summary")
print(summary.round(4).to_string(index=False))

print("\nPairwise class comparisons:")
print(
    tests[[
        "class_a", "class_b", "rank_biserial_a_vs_b",
        "p_value", "BH_q_value"
    ]].round(6).to_string(index=False)
)

print("\nSaved:", output)
print(
    "\nCaution: motif matches use an exploratory threshold. "
    "The 20-shuffle calibration is preliminary and does not establish "
    "TEAD1 occupancy or causal regulation."
)
