import os
import pandas as pd
import matplotlib.pyplot as plt

INPUT = "3_results/tables/motif_prediction_residuals/motif_score_prediction_residual_correlations.csv"
OUTDIR = "3_results/figures/motif_prediction_residuals"

os.makedirs(OUTDIR, exist_ok=True)
df = pd.read_csv(INPUT)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

panels = [
    ("spearman_rho_vs_signed_residual", "BH_q_signed_residual",
     "Motif score vs signed residual",
     "Spearman correlation with signed residual"),
    ("spearman_rho_vs_absolute_error", "BH_q_absolute_error",
     "Motif score vs absolute error",
     "Spearman correlation with absolute error")
]

for ax, (rho_col, q_col, title, xlabel) in zip(axes, panels):
    d = df.sort_values(rho_col).copy()
    colors = ["tab:red" if q < 0.05 else "gray" for q in d[q_col]]
    ax.barh(d["motif"], d[rho_col], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Motif")
    ax.grid(axis="x", linestyle="--", alpha=0.3)

fig.suptitle("Motif scores and AlphaGenome prediction errors", fontsize=13)
fig.tight_layout()

png = os.path.join(OUTDIR, "motif_prediction_residual_correlations.png")
pdf = os.path.join(OUTDIR, "motif_prediction_residual_correlations.pdf")
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

print("Saved:", png)
print("Saved:", pdf)
print("\nRed bars indicate BH-adjusted q < 0.05.")
