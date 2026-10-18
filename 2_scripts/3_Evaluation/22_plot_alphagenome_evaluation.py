from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "3_results" / "tables"
FIGURES = ROOT / "3_results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

evaluation = pd.read_csv(TABLES / "alphagenome_atac_evaluation.csv")
ranking = pd.read_csv(TABLES / "class_ranking_direction.csv")

truth = pd.read_csv(
    ROOT / "1_data_raw/2_benchmark/benchmark_10000_groundtruth.csv"
)
pred = pd.read_csv(
    ROOT / "1_data_raw/3_alphagenome/02_atac_run/"
    "benchmark_10000_alphagenome_atac_pred.csv"
)

keys = ["chr", "start", "end"]
truth["chr"] = truth["chr"].astype(str)
pred["chr"] = pred["chr"].astype(str)

data = truth.merge(
    pred[keys + ["pred_atac_mean"]],
    on=keys,
    how="inner",
    validate="one_to_one",
)

classes = ["RS_gained", "RS_lost", "Stable", "Negative"]
labels = ["RS gained", "RS lost", "Stable", "Negative"]

fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
fig.suptitle(
    "AlphaGenome ATAC-seq Benchmark: WI-38 Replicative Senescence",
    fontsize=15,
    fontweight="bold",
)

# Panel A: correlations with experimental signals
ax = axes[0]
corr = evaluation[evaluation["label"] == "ALL"].set_index("analysis")
metrics = [
    ("pred_vs_young", "Young"),
    ("pred_vs_rs", "Senescent"),
]
pearson = [corr.loc[key, "pearson_r"] for key, _ in metrics]
spearman = [corr.loc[key, "spearman_r"] for key, _ in metrics]
x = np.arange(2)
width = 0.34
ax.bar(x - width / 2, pearson, width, label="Pearson r")
ax.bar(x + width / 2, spearman, width, label="Spearman r")
ax.set_xticks(x, [name for _, name in metrics])
ax.set_ylim(0, 1)
ax.set_ylabel("Correlation")
ax.set_title("A. Agreement with experimental signal")
ax.legend(frameon=False, fontsize=8)
for xpos, value in zip(x - width / 2, pearson):
    ax.text(xpos, value + 0.02, f"{value:.2f}", ha="center", fontsize=8)
for xpos, value in zip(x + width / 2, spearman):
    ax.text(xpos, value + 0.02, f"{value:.2f}", ha="center", fontsize=8)

# Panel B: one-vs-rest ROC-AUC with the better score direction
ax = axes[1]
ranking = ranking.set_index("target_class").loc[classes]
auc = []
for label, row in ranking.iterrows():
    if row["better_direction"] == "higher_score":
        auc.append(row["roc_auc_forward"])
    else:
        auc.append(row["roc_auc_reverse"])
bars = ax.bar(labels, auc)
ax.axhline(0.5, linestyle="--", linewidth=1, label="Chance (0.50)")
ax.set_ylim(0, 1)
ax.set_ylabel("One-vs-rest ROC-AUC")
ax.set_title("B. Class discrimination")
ax.tick_params(axis="x", rotation=25)
ax.legend(frameon=False, fontsize=8)
for bar, value in zip(bars, auc):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        value + 0.02,
        f"{value:.3f}",
        ha="center",
        fontsize=8,
    )

# Panel C: predicted-signal distributions
ax = axes[2]
data["pred_atac_mean"] = pd.to_numeric(
    data["pred_atac_mean"], errors="coerce"
)
box_data = [
    data.loc[data["label"] == label, "pred_atac_mean"].dropna().to_numpy()
    for label in classes
]
ax.boxplot(box_data, tick_labels=labels, showfliers=False)
ax.set_ylabel("Predicted ATAC signal")
ax.set_title("C. Predicted signal by class")
ax.tick_params(axis="x", rotation=25)

fig.text(
    0.5,
    0.015,
    "Benchmark classes are selected from experimental signal extremes or stable peaks; "
    "negative regions are synthetic controls. ROC-AUC is descriptive for this selected benchmark.",
    ha="center",
    fontsize=8,
    wrap=True,
)
fig.tight_layout(rect=[0, 0.07, 1, 0.92])

png = FIGURES / "alphagenome_evaluation_summary.png"
pdf = FIGURES / "alphagenome_evaluation_summary.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

print("Saved:")
print(png)
print(pdf)
print("\nRows joined:", len(data))
print("ROC-AUC values plotted:")
for label, value in zip(labels, auc):
    print(f"  {label}: {value:.4f}")
