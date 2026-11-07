import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

PRED_FILE = "1_data_raw/3_alphagenome/02_atac_run/benchmark_10000_alphagenome_atac_pred.csv"
GT_FILE = "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
MOTIF_FILE = "3_results/tables/unified_motif_analysis/all_nine_motifs_combined.csv"
OUTDIR = "3_results/tables/motif_prediction_residuals"

os.makedirs(OUTDIR, exist_ok=True)

pred = pd.read_csv(PRED_FILE)
gt = pd.read_csv(GT_FILE)
motif = pd.read_csv(MOTIF_FILE)

keys = ["chr", "start", "end", "name", "label"]

# Restrict to the prediction track already present in the prediction file.
if "track_index" in pred.columns:
    tracks = pred["track_index"].unique()
    if len(tracks) != 1:
        raise ValueError(f"Expected one track; found {tracks}")
    print("Using track_index:", tracks[0])

merged = pred.merge(
    gt[keys + ["young_mean", "rs_mean", "log2FC_rs_vs_young"]],
    on=keys,
    how="inner",
    validate="one_to_one"
)

merged["residual_rs"] = merged["pred_atac_mean"] - merged["rs_mean"]
merged["absolute_error_rs"] = merged["residual_rs"].abs()

print("Prediction rows:", len(pred))
print("Ground-truth rows:", len(gt))
print("Matched prediction/ground-truth regions:", len(merged))
print("\nResidual summary:")
print(merged["residual_rs"].describe().to_string())
print("\nMean residual by class:")
print(merged.groupby("label")["residual_rs"].agg(["count", "mean", "median", "std"]).to_string())

# Save region-level residuals.
region_out = os.path.join(OUTDIR, "alphagenome_rs_prediction_residuals.csv")
merged.to_csv(region_out, index=False)
print("\nSaved:", region_out)

# Convert motif table to one row per region, one column per motif.
wide = motif.pivot(index="region_id", columns="motif", values="score_excess")
wide.index.name = "name"
wide = wide.reset_index()

analysis = merged.merge(
    wide,
    on="name",
    how="inner",
    validate="one_to_one"
)

motifs = [c for c in wide.columns if c != "name"]
print("\nRegions with predictions, RS ground truth and motif scores:", len(analysis))
print("Motifs tested:", len(motifs))

rows = []
for m in motifs:
    valid = analysis[[m, "residual_rs", "absolute_error_rs"]].replace(
        [np.inf, -np.inf], np.nan
    ).dropna()

    rho_signed, p_signed = spearmanr(valid[m], valid["residual_rs"])
    rho_abs, p_abs = spearmanr(valid[m], valid["absolute_error_rs"])

    rows.append({
        "motif": m,
        "n_regions": len(valid),
        "spearman_rho_vs_signed_residual": rho_signed,
        "p_vs_signed_residual": p_signed,
        "spearman_rho_vs_absolute_error": rho_abs,
        "p_vs_absolute_error": p_abs
    })

res = pd.DataFrame(rows)
res["BH_q_signed_residual"] = multipletests(
    res["p_vs_signed_residual"], method="fdr_bh"
)[1]
res["BH_q_absolute_error"] = multipletests(
    res["p_vs_absolute_error"], method="fdr_bh"
)[1]

res = res.sort_values("BH_q_signed_residual")
out = os.path.join(OUTDIR, "motif_score_prediction_residual_correlations.csv")
res.to_csv(out, index=False)

print("\nSaved:", out)
print("\nMotif correlations with prediction residuals:")
print(res.to_string(index=False))
print("\nFDR-significant motifs vs signed residual:",
      int((res["BH_q_signed_residual"] < 0.05).sum()))
print("FDR-significant motifs vs absolute error:",
      int((res["BH_q_absolute_error"] < 0.05).sum()))
