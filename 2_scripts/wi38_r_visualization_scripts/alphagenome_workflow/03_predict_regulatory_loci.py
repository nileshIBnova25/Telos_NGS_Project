#!/usr/bin/env python3
"""Run small-scale AlphaGenome reference-sequence predictions around DE genes.

Requires the ontology CURIEs to be chosen after inspecting model metadata:
  export ALPHAGENOME_ONTOLOGY_TERMS='CL:xxxxxxx,UBERON:xxxxxxx'

Runs at most 5 genes per direction by default to keep the first experiment small.
Produces per-gene interactive HTML track plots and a CSV of regional signal summaries.
"""
from pathlib import Path
import os
import re
import sys
import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from alphagenome.data import gene_annotation, genome
from alphagenome.models import dna_client

ROOT = Path.cwd()
OUT_DIR = ROOT / "1_data_raw/4_rnaseq/outputs/alphagenome"
CANDIDATE_PATH = OUT_DIR / "candidate_genes_top20_each_direction.csv"
GTF_URL = "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather"
# First run defaults to the top 5 candidate genes in each direction (up to 10 total).
N_PER_DIRECTION = int(os.environ.get("ALPHAGENOME_N_GENES_PER_DIRECTION", "5"))
SEQUENCE_LENGTH_NAME = os.environ.get("ALPHAGENOME_SEQUENCE_LENGTH", "SEQUENCE_LENGTH_1MB")
OUTPUT_TYPES = [
    dna_client.OutputType.RNA_SEQ,
    dna_client.OutputType.CAGE,
    dna_client.OutputType.DNASE,
    dna_client.OutputType.ATAC,
    dna_client.OutputType.CHIP_HISTONE,
]


def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))


def find_col(df, names):
    norm = {re.sub(r"[^a-z0-9]", "", str(c).lower()): c for c in df.columns}
    for n in names:
        key = re.sub(r"[^a-z0-9]", "", n.lower())
        if key in norm:
            return norm[key]
    return None


def track_rows(track, modality, gene, gene_interval, tss):
    if track is None:
        return [], None
    vals = np.asarray(track.values)
    if vals.ndim == 1:
        vals = vals[:, None]
    meta = track.metadata.copy()
    if vals.shape[1] != len(meta):
        warnings.warn(f"{gene}/{modality}: {vals.shape[1]} data tracks but {len(meta)} metadata rows; skipping summaries")
        return [], None
    # TrackData.interval is genomic coordinate span corresponding to the prediction array.
    tr_interval = getattr(track, "interval", None) or gene_interval
    coord_start = int(tr_interval.start)
    coord_end = int(tr_interval.end)
    n = vals.shape[0]

    def summarize(a, b, j):
        a = max(int(a), coord_start)
        b = min(int(b), coord_end)
        if b <= a:
            return np.nan, np.nan
        i0 = max(0, min(n, int(np.floor((a - coord_start) / (coord_end - coord_start) * n))))
        i1 = max(i0 + 1, min(n, int(np.ceil((b - coord_start) / (coord_end - coord_start) * n))))
        chunk = vals[i0:i1, j]
        if chunk.size == 0:
            return np.nan, np.nan
        return float(np.nanmean(chunk)), float(np.nanmax(chunk))

    rows = []
    # Feature windows: promoter (TSS +/- 2 kb), near-TSS (TSS +/- 50 kb), gene body.
    windows = {
        "promoter_TSS_pm2kb": (tss - 2000, tss + 2000),
        "near_TSS_pm50kb": (tss - 50000, tss + 50000),
        "gene_body": (int(gene_interval.start), int(gene_interval.end)),
    }
    for j in range(vals.shape[1]):
        md = meta.iloc[j].to_dict()
        base = {
            "gene": gene,
            "modality": modality,
            "track_index": j,
            "track_name": md.get("name", ""),
            "biosample_name": md.get("biosample_name", ""),
            "ontology_curie": md.get("ontology_curie", ""),
            "strand": md.get("strand", ""),
            "biosample_type": md.get("biosample_type", ""),
            "data_source": md.get("data_source", ""),
        }
        for wname, (a, b) in windows.items():
            mean, peak = summarize(a, b, j)
            base[wname + "_mean"] = mean
            base[wname + "_max"] = peak
        rows.append(base)
    # Return metadata and arrays for an HTML plot. Downsampling is applied at plot time.
    return rows, (vals, meta, coord_start, coord_end)


def main():
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    ontology_raw = os.environ.get("ALPHAGENOME_ONTOLOGY_TERMS", "").strip()
    if not api_key:
        raise RuntimeError("Set ALPHAGENOME_API_KEY in this terminal; keep it out of scripts and Git.")
    if not ontology_raw:
        raise RuntimeError("Set ALPHAGENOME_ONTOLOGY_TERMS to reviewed ontology CURIEs from the metadata inspection step.")
    terms = [x.strip() for x in ontology_raw.split(",") if x.strip()]
    if not CANDIDATE_PATH.exists():
        raise FileNotFoundError(f"Missing candidate list: {CANDIDATE_PATH}. Run 01_prepare_candidates.py first.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = OUT_DIR / "locus_plots"
    plot_dir.mkdir(exist_ok=True)

    candidates = pd.read_csv(CANDIDATE_PATH)
    gene_col = find_col(candidates, ["gene_for_lookup", "tpm_gene_identifier", "gene_symbol", "gene"])
    direction_col = find_col(candidates, ["direction"])
    if gene_col is None or direction_col is None:
        raise ValueError(f"Candidate file needs gene and direction columns. Found: {list(candidates.columns)}")
    selected = pd.concat([
        candidates[candidates[direction_col].astype(str).eq("Up in RS")].head(N_PER_DIRECTION),
        candidates[candidates[direction_col].astype(str).eq("Down in RS")].head(N_PER_DIRECTION),
    ], ignore_index=True).drop_duplicates(subset=[gene_col])
    if selected.empty:
        raise ValueError("No genes selected from candidate CSV.")

    print("Downloading/loading GENCODE v46 hg38 annotation table (first run requires internet)...")
    # Official AlphaGenome quick-start uses this GENCODE v46 hg38 feather annotation URL.
    gtf = pd.read_feather(GTF_URL)
    model = dna_client.create(api_key)
    seq_len = getattr(dna_client, SEQUENCE_LENGTH_NAME, None)
    if seq_len is None:
        raise ValueError(f"Unknown sequence length {SEQUENCE_LENGTH_NAME!r}. Choose one of {list(dna_client.SUPPORTED_SEQUENCE_LENGTHS)}")

    all_summaries = []
    for gene_number, (_, row) in enumerate(selected.iterrows(), start=1):
        gene = str(row[gene_col]).strip()
        direction = str(row[direction_col])
        if not gene or gene.lower() == "nan":
            continue
        try:
            gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene)
        except Exception as e1:
            try:
                gene_interval = gene_annotation.get_gene_interval(gtf, gene_id=gene)
            except Exception as e2:
                print(f"SKIP {gene}: could not map symbol/ID to hg38 coordinates ({e1}; {e2})")
                continue
        tss = int(gene_interval.start if gene_interval.strand != "-" else gene_interval.end)
        input_interval = gene_interval.resize(seq_len)
        print(f"[{gene_number}/{len(selected)}] {gene} ({direction}): {input_interval.chromosome}:{input_interval.start}-{input_interval.end}; ontology={terms}")
        prediction = model.predict_interval(
            interval=input_interval,
            organism=dna_client.Organism.HOMO_SAPIENS,
            requested_outputs=OUTPUT_TYPES,
            ontology_terms=terms,
        )
        modalities = []
        plot_data = []
        for output_type in OUTPUT_TYPES:
            modality = output_type.name.lower()
            track = getattr(prediction, modality, None)
            if track is None or np.asarray(track.values).size == 0:
                continue
            rows, plotting = track_rows(track, modality, gene, input_interval, tss)
            all_summaries.extend(rows)
            if plotting is not None:
                vals, metadata, start, end = plotting
                modalities.append((modality.upper(), vals, metadata, start, end))
        if not modalities:
            print(f"WARNING: no tracks returned for {gene}. Check selected ontology terms and output metadata.")
            continue

        fig = make_subplots(rows=len(modalities), cols=1, shared_xaxes=True,
                            subplot_titles=[m[0] for m in modalities], vertical_spacing=0.025)
        for r, (modality, vals, metadata, start, end) in enumerate(modalities, start=1):
            step = max(1, int(np.ceil(vals.shape[0] / 12000)))
            x = (start + np.arange(vals.shape[0]) * (end - start) / vals.shape[0] - tss)[::step]
            # Plot up to 8 tracks per modality, using returned metadata to make each trace interpretable.
            ntracks = min(vals.shape[1], 8)
            for j in range(ntracks):
                md = metadata.iloc[j]
                label = str(md.get("biosample_name", md.get("name", f"track_{j}")))
                strand = str(md.get("strand", ""))
                name = f"{label} [{strand}]".strip()
                fig.add_trace(go.Scattergl(
                    x=x, y=vals[::step, j], mode="lines", name=name,
                    legendgroup=f"{modality}_{j}",
                    hovertemplate=f"Gene={gene}<br>Modality={modality}<br>Track={name}<br>Position relative to TSS=%{{x:.0f}} bp<br>Predicted signal=%{{y:.4g}}<extra></extra>",
                ), row=r, col=1)
        fig.add_vline(x=0, line_dash="dash", line_color="black")
        fig.update_layout(
            title=f"AlphaGenome reference-sequence predictions: {gene} ({direction})",
            height=max(600, 260 * len(modalities)),
            hovermode="x unified",
            legend_title="Biosample / strand (up to 8 tracks per modality)",
        )
        fig.update_xaxes(title_text="Position relative to gene TSS (bp)", row=len(modalities), col=1)
        html_path = plot_dir / f"{safe_name(gene)}_{safe_name(direction)}.html"
        fig.write_html(html_path, include_plotlyjs="cdn")
        print(f"  Interactive plot: {html_path}")

    if all_summaries:
        summary = pd.DataFrame(all_summaries)
        summary_path = OUT_DIR / "alphagenome_regulatory_signal_summaries.csv"
        summary.to_csv(summary_path, index=False)
        print(f"\nSaved region-level predicted-signal summaries: {summary_path}")
    else:
        print("No summaries were produced; do not interpret this as absence of regulatory activity. Check gene mapping, selected ontology terms, API output, and metadata.")
    print(f"\nPlots directory: {plot_dir}")
    print("These are reference-sequence model predictions, not measurements from the WI-38 cultures.")

if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
