
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------
# Project paths
# --------------------------------------------------

PROJECT_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..")
)

# Allow input data location to be configured externally
BASE = os.environ.get(
    "TELOS_RNASEQ_BASE",
    "/tmp/hayflick_review"
)

COUNT_FILE = os.path.join(
    BASE, "data", "RNAseq", "count_tables", "RS_raw_counts.txt"
)

META_FILE = os.path.join(
    BASE, "data_processing", "RNAseq", "DEseq", "RS_column_data.txt"
)

OUTDIR = os.path.join(PROJECT_ROOT, "3_results", "figures")
os.makedirs(OUTDIR, exist_ok=True)

if not os.path.isfile(COUNT_FILE):
    raise FileNotFoundError(f"Count file not found: {COUNT_FILE}")

if not os.path.isfile(META_FILE):
    raise FileNotFoundError(f"Metadata file not found: {META_FILE}")

# --------------------------------------------------
# Read counts and metadata
# --------------------------------------------------

counts = pd.read_csv(COUNT_FILE, sep="\t", index_col=0)
meta = pd.read_csv(META_FILE, sep="\t")

counts.columns = counts.columns.str.strip('"')
counts.index = counts.index.astype(str).str.strip('"')
meta.columns = meta.columns.str.strip('"')

for col in ["sample", "TP", "batch"]:
    meta[col] = meta[col].astype(str).str.strip('"')

meta = meta.set_index("sample").loc[counts.columns].reset_index()
meta = meta.rename(columns={"index": "sample"})

# --------------------------------------------------
# Library-size normalisation
# --------------------------------------------------

library_sizes = counts.sum(axis=0)
cpm = counts.div(library_sizes, axis=1) * 1_000_000
log_cpm = np.log2(cpm + 1)

markers = [
    "CDKN1A", "CDKN2A", "MKI67", "PCNA",
    "LMNB1", "IL6", "MMP3", "GLB1"
]

tp_order = [
    "TP1", "TP2", "TP3", "TP4", "TP5",
    "TP6", "TP7", "TP8", "TP9"
]

# --------------------------------------------------
# Plot senescence marker expression
# --------------------------------------------------

fig, axes = plt.subplots(
    2, 4,
    figsize=(16, 8.8),
    sharex=True
)

axes = axes.flatten()

batch_styles = {
    "A": {"color": "tab:blue", "marker": "o"},
    "B": {"color": "tab:orange", "marker": "s"}
}

for ax, gene in zip(axes, markers):

    if gene not in log_cpm.index:
        ax.set_title(gene + " (not found)")
        continue

    values = log_cpm.loc[gene]

    for batch, style in batch_styles.items():

        submeta = meta[meta["batch"] == batch]
        mean_x, mean_y = [], []

        for i, tp in enumerate(tp_order, start=1):

            samples = submeta.loc[
                submeta["TP"] == tp, "sample"
            ].tolist()

            ys = [
                values[s]
                for s in samples
                if s in values.index
            ]

            if not ys:
                continue

            # Replicates are shown individually
            jitter = np.linspace(-0.07, 0.07, len(ys))

            ax.scatter(
                [i + j for j in jitter],
                ys,
                color=style["color"],
                marker=style["marker"],
                s=28,
                alpha=0.8,
                edgecolors="white",
                linewidths=0.35,
                zorder=3
            )

            mean_x.append(i)
            mean_y.append(float(np.mean(ys)))

        # Separate lines per batch
        if mean_x:
            ax.plot(
                mean_x,
                mean_y,
                color=style["color"],
                marker=style["marker"],
                markersize=4,
                linewidth=1.7,
                label=f"Batch {batch} mean",
                zorder=2
            )

    ax.set_title(gene, fontsize=12, pad=8)
    ax.set_ylabel("log2(CPM + 1)", fontsize=9)
    ax.set_xticks(range(1, 10))
    ax.set_xticklabels(tp_order, rotation=0, fontsize=8)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# --------------------------------------------------
# Figure formatting
# --------------------------------------------------

fig.suptitle(
    "Senescence-associated marker expression across the WI-38 time course",
    fontsize=15,
    y=0.99
)

fig.supxlabel("Time point", fontsize=11, y=0.07)

handles = [
    plt.Line2D(
        [0], [0],
        color=style["color"],
        marker=style["marker"],
        linewidth=1.7,
        markersize=5,
        label=f"Batch {batch} mean"
    )
    for batch, style in batch_styles.items()
]

fig.legend(
    handles=handles,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.015),
    ncol=2,
    frameon=False,
    fontsize=9
)

fig.text(
    0.5, 0.045,
    "Points = individual samples; lines = batch-specific means. "
    "TP6 includes batches A (PDL46) and B (PDL45).",
    ha="center",
    fontsize=8
)

fig.subplots_adjust(
    left=0.07,
    right=0.99,
    top=0.87,
    bottom=0.19,
    wspace=0.28,
    hspace=0.38
)

# --------------------------------------------------
# Save figures
# --------------------------------------------------

png = os.path.join(
    OUTDIR, "senescence_marker_timecourse.png"
)

pdf = os.path.join(
    OUTDIR, "senescence_marker_timecourse.pdf"
)

fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")

plt.close(fig)

# --------------------------------------------------
# Summary
# --------------------------------------------------

print("Saved:", png)
print("Saved:", pdf)
print("Samples plotted:", len(meta))
print("Genes plotted:", len(markers))
print("Normalisation: log2(CPM + 1)")
print("Batch A and B means are plotted separately.")
