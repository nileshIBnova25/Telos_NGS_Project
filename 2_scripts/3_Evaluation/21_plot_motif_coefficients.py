from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "3_results/tables/unified_motif_analysis/adjusted_score_excess_tests.csv"
OUTDIR = ROOT / "3_results/figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)

required = {
    "motif", "contrast", "beta_score_excess",
    "ci95_low", "ci95_high", "p_value", "BH_q_value_18_tests"
}
missing = required - set(df.columns)
if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

order = [
    "TEAD1", "CEBPA", "E2F1", "FOS", "JUN",
    "NFKB1", "RELA", "SMAD3", "SMAD4"
]

fig, ax = plt.subplots(figsize=(10, 7))

y_base = np.arange(len(order))
offset = 0.15

for contrast, shift, marker, label in [
    ("RS_gained_vs_Stable", offset, "o", "RS-gained vs Stable"),
    ("RS_lost_vs_Stable", -offset, "s", "RS-lost vs Stable"),
]:
    sub = df[df["contrast"] == contrast].set_index("motif")

    for i, motif in enumerate(order):
        if motif not in sub.index:
            continue

        row = sub.loc[motif]
        beta = float(row["beta_score_excess"])
        low = float(row["ci95_low"])
        high = float(row["ci95_high"])
        q = float(row["BH_q_value_18_tests"])

        # Filled marker denotes BH-adjusted q < 0.05.
        ax.errorbar(
            beta,
            y_base[i] + shift,
            xerr=[[beta - low], [high - beta]],
            fmt=marker,
            markersize=6,
            capsize=3,
            markerfacecolor="black" if q < 0.05 else "white",
            markeredgecolor="black",
            ecolor="black",
            linewidth=1.2,
            label=label if i == 0 else None,
        )

ax.axvline(0, linestyle="--", linewidth=1)
ax.set_yticks(y_base)
ax.set_yticklabels(order)
ax.invert_yaxis()
ax.set_xlabel("Adjusted difference in motif-score excess vs Stable")
ax.set_ylabel("Motif")
ax.set_title("Motif-score excess in replicative-senescence regions")
ax.grid(axis="x", linestyle=":", alpha=0.5)
ax.legend(loc="best", frameon=False)

ax.text(
    0.01, -0.12,
    "Points: adjusted regression coefficients; bars: 95% confidence intervals. "
    "Filled markers: BH q < 0.05.",
    transform=ax.transAxes,
    fontsize=9,
    va="top",
)

fig.tight_layout()

png = OUTDIR / "unified_motif_coefficient_plot.png"
pdf = OUTDIR / "unified_motif_coefficient_plot.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

# Save a concise, report-ready table as well.
table = df.copy()
table["significant_BH_q_lt_005"] = table["BH_q_value_18_tests"] < 0.05
table.to_csv(
    ROOT / "3_results/tables/unified_motif_analysis/"
    "coefficient_plot_results.csv",
    index=False,
)

print("Saved figure:", png)
print("Saved vector PDF:", pdf)
print("Saved coefficient table:",
      ROOT / "3_results/tables/unified_motif_analysis/coefficient_plot_results.csv")
print("\nSignificant comparisons (BH q < 0.05):")
print(
    table.loc[
        table["significant_BH_q_lt_005"],
        ["motif", "contrast", "beta_score_excess",
         "ci95_low", "ci95_high", "BH_q_value_18_tests"]
    ].to_string(index=False)
)
