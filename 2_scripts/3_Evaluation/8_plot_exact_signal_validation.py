import os
import pandas as pd
import matplotlib.pyplot as plt

outdir = (
    "3_results/atac_pipeline/qc/"
    "bigwig_validation/exact_signal_validation"
)

stored_file = (
    "3_results/atac_pipeline/benchmark_10000_replicate_signals.csv"
)
exact_file = os.path.join(outdir, "benchmark_10000_exact_signals.csv")

stored = pd.read_csv(stored_file)
exact = pd.read_csv(exact_file)

tracks = [
    "young_A", "young_B", "young_C",
    "late_TP6_A", "late_TP6_B",
    "late_TP7_A", "late_TP7_B"
]

merged = stored.merge(
    exact,
    on=["chr", "start", "end"],
    suffixes=("_stored", "_exact"),
    validate="one_to_one"
)

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

for i, track in enumerate(tracks):
    ax = axes[i]

    x = merged[f"{track}_stored"]
    y = merged[f"{track}_exact"]

    ax.scatter(x, y, s=5, alpha=0.35)

    low = min(x.min(), y.min())
    high = max(x.max(), y.max())
    ax.plot([low, high], [low, high], linestyle="--")

    ax.set_title(f"{track}\nr = {x.corr(y):.3f}")
    ax.set_xlabel("Stored signal")
    ax.set_ylabel("Exact signal")

axes[-1].axis("off")

fig.suptitle(
    "TELOS: Stored versus Exact ATAC-seq Signal Validation",
    fontsize=14
)
fig.tight_layout()

png_file = os.path.join(outdir, "exact_signal_validation_scatter.png")
pdf_file = os.path.join(outdir, "exact_signal_validation_scatter.pdf")

fig.savefig(png_file, dpi=300, bbox_inches="tight")
fig.savefig(pdf_file, bbox_inches="tight")
plt.close(fig)

print("Saved:", png_file)
print("Saved:", pdf_file)
print("Regions compared per replicate:", len(merged))
print("All seven panels should follow the identity line.")
