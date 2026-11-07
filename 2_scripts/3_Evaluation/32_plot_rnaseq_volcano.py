from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path.cwd()
RNA_FILE = ROOT / "3_results/tables/rnaseq/RS_DESeq2_timecourse_batch_adjusted_unshrunken.csv"
GAINED_FILE = ROOT / "3_results/tables/integration/go_enrichment/RS_gained_genes_with_DESeq2_statistics.csv"
OUT_DIR = ROOT / "3_results/figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

rna = pd.read_csv(RNA_FILE)
gained = pd.read_csv(GAINED_FILE)

rna["minus_log10_padj"] = -np.log10(
    rna["padj"].clip(lower=1e-300)
)

gained_symbols = set(gained["gene_symbol"].dropna())
rna["is_rs_gained"] = rna["gene"].isin(gained_symbols)
rna["significant"] = rna["padj"] < 0.05

# Plot all genes in the background; highlight RS-gained genes.
background = rna[~rna["is_rs_gained"]]
highlight = rna[rna["is_rs_gained"]]

fig, ax = plt.subplots(figsize=(11, 8))

ax.scatter(
    background["log2FoldChange"],
    background["minus_log10_padj"],
    s=9,
    alpha=0.35,
    label="Other genes",
    rasterized=True
)

sig_gained = highlight[highlight["significant"]]
nonsig_gained = highlight[~highlight["significant"]]

ax.scatter(
    nonsig_gained["log2FoldChange"],
    nonsig_gained["minus_log10_padj"],
    s=18,
    alpha=0.8,
    label="RS-gained genes (not significant)",
    rasterized=True
)

ax.scatter(
    sig_gained["log2FoldChange"],
    sig_gained["minus_log10_padj"],
    s=24,
    alpha=0.85,
    label="RS-gained genes (padj < 0.05)",
    rasterized=True
)

ax.axhline(
    -np.log10(0.05),
    linestyle="--",
    linewidth=1,
    label="Adjusted p = 0.05"
)
ax.axvline(0, linestyle=":", linewidth=1)

# Label selected RS-gained genes with the smallest adjusted p-values.
labels = (
    sig_gained.sort_values("padj")
    .drop_duplicates("gene")
    .head(12)
)

for row in labels.itertuples(index=False):
    ax.annotate(
        row.gene,
        (row.log2FoldChange, row.minus_log10_padj),
        xytext=(4, 4),
        textcoords="offset points",
        fontsize=8,
        annotation_clip=True
    )

ax.set_xlabel("DESeq2 time-course log2 fold change", fontsize=11)
ax.set_ylabel("-log10(BH-adjusted p-value)", fontsize=11)
ax.set_title(
    "RNA-seq time-course differential expression\n"
    "RS-gained genes highlighted",
    fontsize=14,
    fontweight="bold"
)
ax.grid(alpha=0.2, linestyle="--")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")

fig.text(
    0.5, 0.015,
    "DESeq2 model: continuous time-course coefficient adjusted for batch. "
    "This is not a direct young-versus-senescent contrast.",
    ha="center",
    fontsize=9
)

fig.tight_layout(rect=[0, 0.045, 1, 1])

png = OUT_DIR / "RNAseq_timecourse_volcano.png"
pdf = OUT_DIR / "RNAseq_timecourse_volcano.pdf"

fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)

print("Genes in DESeq2 table:", len(rna))
print("Genes with padj < 0.05:", int(rna["significant"].sum()))
print("Positive significant genes:", int(
    ((rna["significant"]) & (rna["log2FoldChange"] > 0)).sum()
))
print("Negative significant genes:", int(
    ((rna["significant"]) & (rna["log2FoldChange"] < 0)).sum()
))
print("RS-gained genes highlighted:", len(highlight))
print("RS-gained genes with padj < 0.05:", len(sig_gained))
print("Labels added:", len(labels))
print("Saved:", png)
print("Saved:", pdf)
