#!/usr/bin/env python3
"""Create hover-enabled Plotly HTML versions of the WI-38 RNA-seq plots.

Run from the project root (the directory containing 1_data_raw/):
    python3 1_data_raw/wi38_r_visualization_scripts/05_interactive_html.py

Inputs expected:
  1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx
  1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv
  1_data_raw/4_rnaseq/outputs/visualizations/exploratory_limma_results.csv

Outputs are written to 1_data_raw/4_rnaseq/outputs/visualizations/interactive/.
"""
from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

YOUNG = ["PDL20_TP1_A", "PDL20_TP1_B", "PDL20_TP1_C"]
RS = ["PDL46_TP6_A", "PDL46_TP6_B", "PDL50_TP7_A", "PDL50_TP7_B"]
SAMPLES = YOUNG + RS
COLORS = {"Young": "#f8766d", "RS": "#00bfc4"}
DIRECTION_COLORS = {
    "Not significant": "#a8a8a8",
    "Down in RS": "#4575b4",
    "Up in RS": "#d73027",
}


def first_matching_column(frame: pd.DataFrame, candidates: list[str], what: str) -> str:
    for col in candidates:
        if col in frame.columns:
            return col
    raise ValueError(
        f"Could not find the {what} column. Available columns: {', '.join(map(str, frame.columns))}"
    )


def sample_metadata(sample: str) -> tuple[str, str]:
    condition = "Young" if sample in YOUNG else "RS"
    if sample in YOUNG:
        subgroup = "Young · PDL20"
    elif "PDL46" in sample:
        subgroup = "RS · PDL46 / TP6"
    else:
        subgroup = "RS · PDL50 / TP7"
    return condition, subgroup


def write_plot(fig: go.Figure, path: Path, title: str) -> None:
    fig.update_layout(
        template="plotly_white",
        title=dict(text=title, x=0.02, xanchor="left"),
        hoverlabel=dict(namelength=-1),
        margin=dict(l=80, r=35, t=100, b=75),
        font=dict(family="Arial, sans-serif", size=13),
    )
    fig.write_html(
        path,
        include_plotlyjs="directory",
        full_html=True,
        config={"responsive": True, "displaylogo": False, "scrollZoom": True},
    )
    print(f"Created: {path}")


def make_volcano(results_path: Path, out_dir: Path) -> dict:
    res = pd.read_csv(results_path)
    gene_col = first_matching_column(res, ["gene", "gene_id_or_symbol", "gene_symbol"], "gene")
    fc_col = first_matching_column(res, ["logFC", "log2FC", "log2FoldChange"], "log fold change")
    p_col = first_matching_column(res, ["P.Value", "p_value", "pvalue", "p_val"], "nominal p-value")
    padj_col = first_matching_column(res, ["adj.P.Val", "padj", "FDR", "q_value"], "adjusted p-value")
    if "direction" not in res:
        fc = pd.to_numeric(res[fc_col], errors="coerce")
        padj = pd.to_numeric(res[padj_col], errors="coerce")
        res["direction"] = np.where(
            (padj < 0.05) & (fc >= 1), "Up in RS",
            np.where((padj < 0.05) & (fc <= -1), "Down in RS", "Not significant")
        )
    if "neglog10_p" in res:
        res["_y"] = pd.to_numeric(res["neglog10_p"], errors="coerce")
    else:
        p = pd.to_numeric(res[p_col], errors="coerce").clip(lower=1e-300)
        res["_y"] = -np.log10(p)
    for col in [fc_col, p_col, padj_col]:
        res[col] = pd.to_numeric(res[col], errors="coerce")
    res = res.replace([np.inf, -np.inf], np.nan).dropna(subset=[gene_col, fc_col, "_y", p_col, padj_col])
    fig = go.Figure()
    hover = np.column_stack([
        res[gene_col].astype(str), res[p_col].to_numpy(), res[padj_col].to_numpy(),
        res.get("AveExpr", pd.Series([np.nan] * len(res))).to_numpy(),
        res.get("t", pd.Series([np.nan] * len(res))).to_numpy(),
    ])
    # Plot non-significant first so significant points sit above them visually.
    for direction in ["Not significant", "Down in RS", "Up in RS"]:
        sub = res[res["direction"] == direction]
        if sub.empty:
            continue
        custom = np.column_stack([
            sub[gene_col].astype(str).to_numpy(),
            sub[p_col].to_numpy(), sub[padj_col].to_numpy(),
            sub.get("AveExpr", pd.Series([np.nan] * len(sub))).to_numpy(),
            sub.get("t", pd.Series([np.nan] * len(sub))).to_numpy(),
        ])
        fig.add_trace(go.Scattergl(
            x=sub[fc_col], y=sub["_y"], mode="markers", name=direction,
            marker=dict(color=DIRECTION_COLORS.get(direction, "#888888"), size=5, opacity=0.72),
            customdata=custom,
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>Class: " + html.escape(direction) +
                "<br>log2 fold change (RS vs Young): %{x:.3f}" +
                "<br>−log10 nominal p-value: %{y:.3f}" +
                "<br>Nominal p-value: %{customdata[1]:.3g}" +
                "<br>Adjusted p-value: %{customdata[2]:.3g}" +
                "<br>Average log-expression (AveExpr): %{customdata[3]:.3f}" +
                "<br>Moderated t: %{customdata[4]:.3f}<extra></extra>"
            ),
        ))
    fig.add_vline(x=-1, line_dash="dash", line_color="#444", line_width=1)
    fig.add_vline(x=1, line_dash="dash", line_color="#444", line_width=1)
    fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="#444", line_width=1)
    fig.update_layout(
        xaxis_title="log2 fold change (RS vs Young)",
        yaxis_title="−log10 nominal p-value",
        legend_title="Exploratory classification",
        hovermode="closest",
    )
    write_plot(fig, out_dir / "volcano_interactive.html",
               "WI-38 young vs replicative senescence — interactive volcano plot")
    fc = res[fc_col]
    padj = res[padj_col]
    up_rows = res[res["direction"] == "Up in RS"].sort_values(padj_col)
    down_rows = res[res["direction"] == "Down in RS"].sort_values(padj_col)
    return {
        "n_genes": int(len(res)),
        "n_up": int(((padj < 0.05) & (fc >= 1)).sum()),
        "n_down": int(((padj < 0.05) & (fc <= -1)).sum()),
        "n_fdr05": int((padj < 0.05).sum()),
        "labels": res.sort_values(padj_col)[gene_col].astype(str).head(12).tolist(),
        "top_up": up_rows[gene_col].astype(str).head(8).tolist(),
        "top_down": down_rows[gene_col].astype(str).head(8).tolist(),
    }


def make_heatmap(gt_path: Path, tpm: pd.DataFrame, out_dir: Path) -> dict:
    gt = pd.read_csv(gt_path, low_memory=False)
    gene_col = first_matching_column(
        gt, ["gene_id_or_symbol", "gene", "gene_symbol", "gene_name", "Gene", "GeneSymbol", "gene_id", "Gene ID"], "ground-truth gene ID"
    )
    fc_col = first_matching_column(
        gt, ["log2FC_rs_vs_young", "log2FC", "log2FoldChange", "logFC", "log2_fc_rs_vs_young"], "ground-truth log2 fold change"
    )
    if "gene" not in tpm.columns:
        raise ValueError(f"TPM sheet must contain a 'gene' column. Available: {', '.join(map(str, tpm.columns))}")
    missing = [s for s in SAMPLES if s not in tpm.columns]
    if missing:
        raise ValueError(f"TPM sheet is missing expected sample columns: {missing}")
    gt[gene_col] = gt[gene_col].astype("string").str.strip()
    gt[fc_col] = pd.to_numeric(gt[fc_col], errors="coerce")
    gt = gt.dropna(subset=[gene_col, fc_col])
    gt = gt[gt[gene_col].astype(str).ne("")].drop_duplicates(gene_col)
    gt = gt.assign(_absfc=gt[fc_col].abs()).sort_values("_absfc", ascending=False).head(100)

    expr = tpm.loc[:, ["gene"] + SAMPLES].copy()
    expr["gene"] = expr["gene"].astype("string").str.strip()
    expr = expr.dropna(subset=["gene"]).drop_duplicates("gene", keep="first").set_index("gene")
    gt_ids = gt[gene_col].astype(str).tolist()
    matched_ids = [gene for gene in gt_ids if gene in expr.index]
    if len(matched_ids) < 2:
        raise ValueError(
            f"Only {len(matched_ids)} top ground-truth genes matched the TPM gene column. "
            "Check whether the two files use the same gene identifiers."
        )
    fc_map = gt.set_index(gt[gene_col].astype(str))[fc_col].to_dict()
    matched_expr = expr.loc[matched_ids, SAMPLES].apply(pd.to_numeric, errors="coerce")
    matched_expr = matched_expr.replace([np.inf, -np.inf], np.nan)
    matched_expr = matched_expr.dropna(axis=0, how="any")
    log_expr = np.log2(matched_expr.clip(lower=0) + 0.001)
    row_sd = log_expr.std(axis=1, ddof=1).replace(0, np.nan)
    z = log_expr.sub(log_expr.mean(axis=1), axis=0).div(row_sd, axis=0).replace([np.inf, -np.inf], np.nan).dropna(axis=0)
    if z.empty:
        raise ValueError("No usable genes remained for the heatmap after checking variance and missing values.")
    # Add exact underlying measurements to every heatmap cell's hover payload.
    custom = np.empty((z.shape[0], z.shape[1], 5), dtype=object)
    for i, gene in enumerate(z.index.astype(str)):
        for j, sample in enumerate(z.columns):
            condition, subgroup = sample_metadata(sample)
            custom[i, j, :] = [gene, sample, condition, matched_expr.loc[gene, sample], fc_map.get(gene, np.nan)]
    fig = go.Figure(go.Heatmap(
        z=z.to_numpy(), x=z.columns.astype(str).tolist(), y=z.index.astype(str).tolist(),
        customdata=custom, colorscale="RdYlBu_r", zmid=0,
        colorbar=dict(title="Row Z-score", thickness=16),
        hovertemplate=(
            "<b>Gene: %{customdata[0]}</b><br>Sample: %{customdata[1]}" +
            "<br>Condition: %{customdata[2]}" +
            "<br>TPM: %{customdata[3]:.4g}" +
            "<br>Ground-truth log2FC (RS vs Young): %{customdata[4]:.3f}" +
            "<br>Row Z-score of log2(TPM + 0.001): %{z:.3f}<extra></extra>"
        ),
        xgap=1, ygap=1,
    ))
    fig.update_layout(
        xaxis_title="Sample", yaxis_title="Gene (top 100 by absolute ground-truth log2FC)",
        height=max(850, min(1900, 350 + 17 * len(z))),
        margin=dict(l=150, r=45, t=100, b=100),
    )
    fig.update_yaxes(autorange="reversed", tickfont=dict(size=9))
    fig.update_xaxes(side="bottom", tickangle=-30)
    write_plot(fig, out_dir / "heatmap_interactive.html",
               "Top changing WI-38 genes — hover a cell for gene, sample, TPM and Z-score")
    return {"n_selected": int(len(gt_ids)), "n_matched": int(len(matched_ids)), "n_plotted": int(len(z))}


def make_pca(tpm: pd.DataFrame, out_dir: Path) -> dict:
    if "gene" not in tpm.columns:
        raise ValueError("TPM sheet must contain a 'gene' column.")
    missing = [s for s in SAMPLES if s not in tpm.columns]
    if missing:
        raise ValueError(f"TPM sheet is missing expected sample columns: {missing}")
    expr = tpm.loc[:, ["gene"] + SAMPLES].copy()
    expr["gene"] = expr["gene"].astype("string").str.strip()
    expr = expr.dropna(subset=["gene"]).drop_duplicates("gene", keep="first")
    mat = expr[SAMPLES].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    good_rows = np.isfinite(mat).all(axis=1) & (np.nanmean(mat, axis=1) > 1)
    mat = mat[good_rows, :]
    genes = expr.loc[good_rows, "gene"].astype(str).to_numpy()
    if mat.shape[0] < 3:
        raise ValueError(f"Too few genes remain for PCA after mean TPM > 1 filtering: {mat.shape[0]}")
    logmat = np.log2(np.clip(mat, 0, None) + 0.001)
    variances = np.var(logmat, axis=1, ddof=1)
    valid = np.isfinite(variances)
    logmat, variances = logmat[valid], variances[valid]
    genes = genes[valid]
    top_n = min(1000, len(variances))
    selected = np.argsort(variances)[::-1][:top_n]
    sample_matrix = logmat[selected, :].T  # samples x genes
    centered = sample_matrix - sample_matrix.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(centered, full_matrices=False)
    scores = u * singular
    pc_var = singular**2 / max(1, sample_matrix.shape[0] - 1)
    ratio = pc_var / pc_var.sum() * 100
    fig = go.Figure()
    for subgroup in ["Young · PDL20", "RS · PDL46 / TP6", "RS · PDL50 / TP7"]:
        indices = []
        for i, sample in enumerate(SAMPLES):
            _, group = sample_metadata(sample)
            if group == subgroup:
                indices.append(i)
        if not indices:
            continue
        custom = []
        for i in indices:
            sample = SAMPLES[i]
            condition, sub = sample_metadata(sample)
            custom.append([sample, condition, sub, scores[i, 0], scores[i, 1]])
        condition = "Young" if subgroup.startswith("Young") else "RS"
        symbol = "circle" if condition == "Young" else ("triangle-up" if "PDL46" in subgroup else "diamond")
        fig.add_trace(go.Scatter(
            x=scores[indices, 0], y=scores[indices, 1], mode="markers",
            name=subgroup, marker=dict(size=13, color=COLORS[condition], symbol=symbol,
                                       line=dict(color="#333333", width=0.7)),
            customdata=np.asarray(custom, dtype=object),
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>Condition: %{customdata[1]}" +
                "<br>Sample group: %{customdata[2]}" +
                "<br>PC1 score: %{customdata[3]:.3f}" +
                "<br>PC2 score: %{customdata[4]:.3f}<extra></extra>"
            ),
        ))
    fig.update_layout(
        xaxis_title=f"PC1 ({ratio[0]:.1f}% variance)",
        yaxis_title=f"PC2 ({ratio[1]:.1f}% variance)",
        legend_title="Sample group", hovermode="closest",
        annotations=[dict(
            x=0.01, y=-0.16, xref="paper", yref="paper", showarrow=False,
            text=f"PCA of top {top_n} most variable genes using log2(TPM + 0.001). Hover points for sample IDs.",
            align="left", font=dict(size=11, color="#555555")
        )],
    )
    write_plot(fig, out_dir / "pca_interactive.html", "PCA of WI-38 RNA expression — interactive sample explorer")
    return {"n_genes": int(top_n), "pc1": float(ratio[0]), "pc2": float(ratio[1]), "total_pc12": float(ratio[:2].sum())}


def write_index(out_dir: Path, volcano: dict, heatmap: dict, pca: dict) -> None:
    # Findings below are deliberately framed as exploratory interpretation, not proof.
    top_genes = ", ".join(html.escape(g) for g in volcano["labels"][:12])
    up_genes = ", ".join(html.escape(g) for g in volcano["top_up"]) or "No genes met the plotted up-in-RS rule."
    down_genes = ", ".join(html.escape(g) for g in volcano["top_down"]) or "No genes met the plotted down-in-RS rule."
    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>WI-38 interactive RNA-seq explorer</title>
<style>
body{{font:16px/1.55 system-ui,-apple-system,Segoe UI,sans-serif;max-width:1050px;margin:40px auto;padding:0 22px;color:#202124}}
h1{{line-height:1.2}} h2{{margin-top:1.7em}} .links a{{display:inline-block;margin:8px 10px 8px 0;padding:10px 13px;border:1px solid #c9cdd2;border-radius:8px;color:#174ea6;text-decoration:none}}
.note{{border-left:4px solid #777;padding:8px 16px;background:#f7f7f7}} code{{background:#f2f2f2;padding:2px 4px}}
</style></head><body>
<h1>WI-38 RNA-seq: interactive plot explorer</h1>
<p>Open a plot, then hover over points or heatmap cells for gene/sample details. Use the Plotly toolbar to zoom, pan, reset axes, or save a plot as a PNG.</p>
<div class="links"><a href="volcano_interactive.html">Open interactive volcano plot</a><a href="heatmap_interactive.html">Open interactive heatmap</a><a href="pca_interactive.html">Open interactive PCA</a></div>
<h2>What these plots suggest</h2>
<h3>Sample structure (PCA)</h3>
<p>The submitted static PCA shows the three young PDL20 samples grouped closely together, while the four late/RS samples are displaced from them. This is consistent with a substantial transcriptome difference between young and late samples. The RS samples also split into a PDL46/TP6 pair and a PDL50/TP7 pair, suggesting appreciable heterogeneity within the late samples. Because passage/timepoint and condition may be entangled, this split should not automatically be called a biological subtype or a technical batch effect. PC1 explained about 72.9% and PC2 about 11.3% in the original static plot; regenerated values are shown in the interactive PCA axes.</p>
<h3>Expression heatmap</h3>
<p>The heatmap displays row Z-scores, not absolute expression: colors indicate each gene's relative expression across the seven samples. The broad young-versus-late pattern supports group-level differences among the selected genes, while sample-specific bright cells indicate genes that vary particularly strongly in one replicate. Hover a cell to see its gene, sample, TPM, ground-truth log2 fold change, and row Z-score. The top {heatmap['n_plotted']} usable rows were plotted; {heatmap['n_matched']} of the {heatmap['n_selected']} selected ground-truth genes matched to the TPM table.</p>
<h3>Volcano plot</h3>
<p>The volcano plot uses the previously generated exploratory limma results from log2(TPM + 0.001). In that model, positive log2 fold change indicates higher expression in RS and negative log2 fold change indicates lower expression in RS. Top genes classified as <strong>Up in RS</strong> by adjusted p-value &lt; 0.05 and log2FC ≥ 1: {up_genes}. Top genes classified as <strong>Down in RS</strong> by adjusted p-value &lt; 0.05 and log2FC ≤ −1: {down_genes}. The adjusted-p-value-ranked genes in the input results include: {top_genes}. Inspect all genes through hover rather than relying only on these top names. Very large fold changes can arise when expression is near zero in one group, so check per-sample TPM before biological interpretation.</p>
<p>In the generated table, {volcano['n_fdr05']} genes have adjusted p-value &lt; 0.05, including {volcano['n_up']} meeting the additional exploratory threshold log2FC ≥ 1 and {volcano['n_down']} meeting log2FC ≤ −1. These are summaries of the TPM-based limma model and are not definitive count-based differential-expression results.</p>
<h2>Scientific caveat</h2>
<div class="note"><strong>Exploratory only.</strong> TPM values are not raw RNA-seq counts. The p-values in the volcano plot should not be presented as definitive differential-expression evidence. For formal inference, use the appropriate raw counts or estimated counts with a method such as DESeq2, edgeR, or limma-voom, and verify the sample design (including passage/timepoint and any batch variables).</div>
<h2>Files</h2><p>All HTML files in this folder are interactive. Keep the accompanying <code>plotly.min.js</code> file beside the HTML pages so they can work offline.</p>
</body></html>"""
    (out_dir / "index.html").write_text(page, encoding="utf-8")
    # A short machine-readable summary for record-keeping.
    summary = pd.DataFrame([
        {"plot": "volcano", "metric": "genes in result", "value": volcano["n_genes"]},
        {"plot": "volcano", "metric": "adjusted p < 0.05", "value": volcano["n_fdr05"]},
        {"plot": "volcano", "metric": "up in RS (adj p < 0.05 and log2FC >= 1)", "value": volcano["n_up"]},
        {"plot": "volcano", "metric": "down in RS (adj p < 0.05 and log2FC <= -1)", "value": volcano["n_down"]},
        {"plot": "heatmap", "metric": "matched selected genes", "value": heatmap["n_matched"]},
        {"plot": "heatmap", "metric": "genes plotted", "value": heatmap["n_plotted"]},
        {"plot": "PCA", "metric": "PC1 variance percent", "value": round(pca["pc1"], 2)},
        {"plot": "PCA", "metric": "PC2 variance percent", "value": round(pca["pc2"], 2)},
    ])
    summary.to_csv(out_dir / "interactive_summary.csv", index=False)
    print(f"Created: {out_dir / 'index.html'}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd(), help="Project root containing 1_data_raw/ (default: current directory)")
    args = parser.parse_args()
    root = args.project_root.expanduser().resolve()
    gt_path = root / "1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv"
    tpm_path = root / "1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx"
    vis_dir = root / "1_data_raw/4_rnaseq/outputs/visualizations"
    results_path = vis_dir / "exploratory_limma_results.csv"
    out_dir = vis_dir / "interactive"
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in [gt_path, tpm_path, results_path]:
        if not path.is_file():
            print(f"ERROR: required input not found: {path}", file=sys.stderr)
            if path == results_path:
                print("Run ./1_data_raw/wi38_r_visualization_scripts/run_all.sh from the project root first.", file=sys.stderr)
            return 2
    try:
        tpm = pd.read_excel(tpm_path, sheet_name="RS_TC_raw_TPM", engine="calamine")
        volcano_info = make_volcano(results_path, out_dir)
        heatmap_info = make_heatmap(gt_path, tpm, out_dir)
        pca_info = make_pca(tpm, out_dir)
        write_index(out_dir, volcano_info, heatmap_info, pca_info)
    except Exception as exc:
        print(f"ERROR while creating interactive plots: {exc}", file=sys.stderr)
        return 1
    print("\nDone. Open this file in a browser:")
    print(f"  {out_dir / 'index.html'}")
    print("These are self-contained HTML plot pages using a shared local Plotly JavaScript file.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
