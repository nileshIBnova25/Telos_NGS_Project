from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec

ROOT = Path.cwd()
BASE = ROOT / "3_results/tables/integration/go_enrichment"
SENS = BASE / "no_protocadherins"
OUT = ROOT / "3_results/figures"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "Original": BASE / "GO_BP_enrichment_all_tested_terms.csv",
    "Without protocadherins":
        SENS / "GO_BP_enrichment_all_tested_terms.csv",
}

PANELS = [
    ("Original", "RS_gained_all",
     "All RS-gained genes", "Original analysis"),
    ("Original", "RS_gained_concordant",
     "Concordant RS-gained genes", "Original analysis"),
    ("Without protocadherins", "RS_gained_all",
     "All RS-gained genes", "Excluding protocadherins"),
    ("Without protocadherins", "RS_gained_concordant",
     "Concordant RS-gained genes", "Excluding protocadherins"),
]

# Load significant GO terms for each panel.
panel_data = []

for source, analysis, label, subtitle in PANELS:
    df = pd.read_csv(FILES[source])

    df = df[
        (df["analysis"] == analysis)
        & (df["adjusted_p_value"] < 0.05)
        & (df["overlap_n"] >= 3)
    ].copy()

    df["minus_log10_q"] = -np.log10(
        df["adjusted_p_value"].clip(lower=1e-300)
    )

    df = df.sort_values("fold_enrichment").reset_index(drop=True)
    panel_data.append(df)

nonempty = [df for df in panel_data if not df.empty]

if not nonempty:
    raise ValueError("No significant GO terms found in the input files.")

all_data = pd.concat(nonempty, ignore_index=True)

norm = Normalize(
    vmin=all_data["minus_log10_q"].min(),
    vmax=all_data["minus_log10_q"].max()
)

# Figure with a separate label column and plotting column per panel.
fig = plt.figure(figsize=(19, 14))
outer = GridSpec(
    2, 2,
    figure=fig,
    left=0.055,
    right=0.83,
    bottom=0.12,
    top=0.88,
    wspace=0.18,
    hspace=0.38
)

scatter_axes = []

for panel_index, (panel, df) in enumerate(zip(PANELS, panel_data)):
    source, analysis, label, subtitle = panel
    row = panel_index // 2
    col = panel_index % 2

    inner = outer[row, col].subgridspec(
        1, 2,
        width_ratios=[2.25, 1.45],
        wspace=0.03
    )

    label_ax = fig.add_subplot(inner[0, 0])
    ax = fig.add_subplot(inner[0, 1])
    scatter_axes.append(ax)

    label_ax.set_title(
        f"{label}\n{subtitle}",
        fontsize=12.5,
        fontweight="bold",
        pad=12
    )

    if df.empty:
        label_ax.text(
            0.95, 0.5, "No significant terms",
            ha="right", va="center",
            transform=label_ax.transAxes,
            fontsize=10
        )
        label_ax.axis("off")
        ax.axis("off")
        continue

    n = len(df)
    y = np.arange(n)

    # Dedicated text axis: labels never extend into the scatter plot.
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(-0.5, n - 0.5)
    label_ax.set_yticks([])
    label_ax.set_xticks([])

    for spine in label_ax.spines.values():
        spine.set_visible(False)

    for i, term in enumerate(df["term"]):
        label_ax.text(
            0.98, i, term,
            ha="right",
            va="center",
            fontsize=9.5,
            wrap=False,
            clip_on=True
        )

    ax.scatter(
        df["fold_enrichment"],
        y,
        c=df["minus_log10_q"],
        cmap="viridis",
        norm=norm,
        s=35 + 8 * df["overlap_n"],
        edgecolors="black",
        linewidths=0.45,
        alpha=0.95
    )

    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Fold enrichment", fontsize=10, labelpad=7)
    ax.tick_params(axis="x", labelsize=9)
    ax.grid(axis="x", linestyle="--", linewidth=0.6, alpha=0.3)
    ax.set_axisbelow(True)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    xmax = df["fold_enrichment"].max()
    ax.set_xlim(0, xmax * 1.18 if xmax > 0 else 1)

# Shared significance colour scale in a separate figure margin.
sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
sm.set_array([])

cax = fig.add_axes([0.875, 0.32, 0.017, 0.40])
cbar = fig.colorbar(sm, cax=cax)
cbar.set_label("-log10(BH-adjusted q-value)", fontsize=10, labelpad=9)
cbar.ax.tick_params(labelsize=9)

# Bubble-size legend.
sizes = sorted(set([
    int(all_data["overlap_n"].min()),
    int(all_data["overlap_n"].median()),
    int(all_data["overlap_n"].max())
]))

handles = [
    plt.scatter(
        [], [],
        s=35 + 8 * size,
        facecolor="lightgray",
        edgecolor="black",
        linewidth=0.5,
        label=str(size)
    )
    for size in sizes
]

fig.legend(
    handles=handles,
    title="Overlapping genes",
    loc="center left",
    bbox_to_anchor=(0.855, 0.20),
    frameon=False,
    fontsize=9,
    title_fontsize=9
)

fig.suptitle(
    "GO Biological Process enrichment in RS-gained genes",
    fontsize=17,
    fontweight="bold",
    y=0.96
)

fig.text(
    0.5, 0.045,
    "Only terms with BH-adjusted q < 0.05 and at least 3 overlapping genes "
    "are shown. Dot size represents overlap count; colour represents "
    "statistical significance. Lower panels exclude PCDH-family genes "
    "from both the selected sets and background.",
    ha="center",
    va="center",
    fontsize=9,
    wrap=True
)

png = OUT / "GO_BP_enrichment_comparison.png"
pdf = OUT / "GO_BP_enrichment_comparison.pdf"

fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("Saved PNG:", png)
print("Saved PDF:", pdf)
print("\nSignificant terms plotted:")
for panel, df in zip(PANELS, panel_data):
    print(f"{panel[2]} | {panel[3]}: {len(df)} terms")
