#!/usr/bin/env bash
# Build a curated, presentation-friendly ZIP of the important Telos_NGS_Project
# result tables, figures, enrichment outputs, reports, and analysis scripts.
# Raw sequencing data, genome indexes, large tracks, and .git are deliberately excluded.
set -Eeuo pipefail

ROOT="${HOME}/Sem3/BDBP307_Omics_Lab/ProjectNGS/Telos_NGS_Project"
if [[ ! -d "$ROOT" ]]; then
  echo "ERROR: Project directory not found: $ROOT" >&2
  echo "Edit ROOT near the top of this script if your project is stored elsewhere." >&2
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
PACKAGE_NAME="Telos_NGS_Presentation_Assets_${STAMP}"
REPORTS_DIR="$ROOT/reports"
PKG="$REPORTS_DIR/$PACKAGE_NAME"
ZIP_PATH="$REPORTS_DIR/${PACKAGE_NAME}.zip"
RNA="$ROOT/1_data_raw/4_rnaseq"
OUT="$RNA/outputs"

mkdir -p \
  "$PKG/00_Project_Docs" \
  "$PKG/01_RNAseq/00_Presentation_Figures" \
  "$PKG/01_RNAseq/01_QC_and_Input_Inspection" \
  "$PKG/01_RNAseq/02_Exploratory_Visualizations" \
  "$PKG/01_RNAseq/03_Expression_Tables" \
  "$PKG/01_RNAseq/04_Downstream_Package" \
  "$PKG/02_YAP_TEAD_LISA/01_Gene_Sets_and_Analysis" \
  "$PKG/03_Functional_Enrichment" \
  "$PKG/04_AlphaGenome" \
  "$PKG/05_ATAC_and_Multiomics/Discovered_Project_Files" \
  "$PKG/06_Methods_and_Scripts" \
  "$PKG/07_Other_Downloaded_Results"

# Result/document types safe to bundle. Large raw instrument/genome files are excluded.
allowed_ext() {
  case "${1##*.}" in
    csv|tsv|txt|pdf|png|jpg|jpeg|svg|html|md|docx|xlsx|xml|json|R|r|py|sh|bed|narrowPeak|gmt|yaml|yml) return 0 ;;
    *) return 1 ;;
  esac
}

copy_file() {
  local src="$1" dst="$2"
  [[ -f "$src" ]] || return 0
  mkdir -p "$(dirname "$dst")"
  cp -p "$src" "$dst"
}

copy_tree_results() {
  local src="$1" dst="$2" f rel
  [[ -d "$src" ]] || return 0
  while IFS= read -r -d '' f; do
    allowed_ext "$f" || continue
    rel="${f#"$src"/}"
    copy_file "$f" "$dst/$rel"
  done < <(find "$src" -type f -print0)
}

# Copy the known RNA-seq deliverable directories, keeping their original filenames.
copy_tree_results "$OUT/rnaseq_downstream_package" "$PKG/01_RNAseq/04_Downstream_Package"
copy_tree_results "$OUT/visualizations" "$PKG/01_RNAseq/02_Exploratory_Visualizations"
copy_tree_results "$OUT/yap_tead_rna" "$PKG/02_YAP_TEAD_LISA/01_Gene_Sets_and_Analysis"
copy_tree_results "$OUT/alphagenome" "$PKG/04_AlphaGenome"
copy_tree_results "$RNA/qc" "$PKG/01_RNAseq/01_QC_and_Input_Inspection/QC"
copy_tree_results "$RNA/inspection" "$PKG/01_RNAseq/01_QC_and_Input_Inspection/Input_Inspection"

# Any direct, output-root result tables (e.g. ground-truth TPM / integrated CSVs).
if [[ -d "$OUT" ]]; then
  while IFS= read -r -d '' f; do
    allowed_ext "$f" || continue
    copy_file "$f" "$PKG/01_RNAseq/03_Expression_Tables/$(basename "$f")"
  done < <(find "$OUT" -maxdepth 1 -type f -print0)
fi

# Copy root-level project documentation and report files (not raw data).
for f in "$ROOT"/README* "$ROOT"/*.md "$ROOT"/*.txt "$REPORTS_DIR"/*.docx "$REPORTS_DIR"/*.pdf "$REPORTS_DIR"/*.md "$REPORTS_DIR"/*.txt; do
  [[ -f "$f" ]] || continue
  copy_file "$f" "$PKG/00_Project_Docs/$(basename "$f")"
done

# Include analysis scripts and method READMEs, preserving folder structure.
if [[ -d "$ROOT/2_scripts" ]]; then
  while IFS= read -r -d '' f; do
    case "$f" in
      *.R|*.r|*.py|*.sh|*.md|*.txt) ;;
      *) continue ;;
    esac
    rel="${f#"$ROOT/2_scripts/"}"
    copy_file "$f" "$PKG/06_Methods_and_Scripts/$rel"
  done < <(find "$ROOT/2_scripts" -type f -print0)
fi

# Make presentation-friendly copies of the main new RNA figures with clearer names.
copy_file "$OUT/rnaseq_downstream_package/figure_01_GSEA_Hallmark_directional.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig01_Hallmark_GSEA_Directional.png"
copy_file "$OUT/rnaseq_downstream_package/figure_01_GSEA_Hallmark_directional.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig01_Hallmark_GSEA_Directional.pdf"
copy_file "$OUT/rnaseq_downstream_package/figure_02_candidate_gene_trajectories.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig02_Candidate_Gene_Trajectories.png"
copy_file "$OUT/rnaseq_downstream_package/figure_02_candidate_gene_trajectories.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig02_Candidate_Gene_Trajectories.pdf"
copy_file "$OUT/rnaseq_downstream_package/figure_03_endpoint_vs_time_trend.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig03_Endpoint_TPM_vs_Time_Trend.png"
copy_file "$OUT/rnaseq_downstream_package/figure_03_endpoint_vs_time_trend.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig03_Endpoint_TPM_vs_Time_Trend.pdf"
copy_file "$OUT/rnaseq_downstream_package/figure_04_alphagenome_predRNA_vs_observed_TPM.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig04_AlphaGenome_Predicted_RNA_vs_TPM.png"
copy_file "$OUT/rnaseq_downstream_package/figure_04_alphagenome_predRNA_vs_observed_TPM.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig04_AlphaGenome_Predicted_RNA_vs_TPM.pdf"
copy_file "$OUT/visualizations/volcano_plot.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig05_Exploratory_Volcano.png"
copy_file "$OUT/visualizations/volcano_plot.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig05_Exploratory_Volcano.pdf"
copy_file "$OUT/visualizations/heatmap_top100_genes.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig06_Heatmap_Top100_Genes.png"
copy_file "$OUT/visualizations/heatmap_top100_genes.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig06_Heatmap_Top100_Genes.pdf"
copy_file "$OUT/visualizations/pca_young_vs_rs.png" "$PKG/01_RNAseq/00_Presentation_Figures/Fig07_PCA_Young_vs_Late.png"
copy_file "$OUT/visualizations/pca_young_vs_rs.pdf" "$PKG/01_RNAseq/00_Presentation_Figures/Fig07_PCA_Young_vs_Late.pdf"

# Friendly filenames for the main report/PPT tables and gene lists.
copy_file "$OUT/rnaseq_downstream_package/GSEA_Hallmark_top_directional_terms.csv" "$PKG/03_Functional_Enrichment/Table01_Hallmark_GSEA_Top_Directional_Terms.csv"
copy_file "$OUT/rnaseq_downstream_package/candidate_gene_panel_summary.csv" "$PKG/01_RNAseq/03_Expression_Tables/Table02_Candidate_Regulators_and_Senescence_Genes.csv"
copy_file "$OUT/rnaseq_downstream_package/time_trend_vs_endpoint_agreement_metrics.csv" "$PKG/01_RNAseq/03_Expression_Tables/Table03_Time_Trend_vs_Endpoint_Agreement.csv"
copy_file "$OUT/rnaseq_downstream_package/alphagenome_predRNA_vs_observed_TPM_metrics.csv" "$PKG/04_AlphaGenome/Table04_Predicted_RNA_Signal_vs_Observed_TPM_Metrics.csv"
copy_file "$OUT/rnaseq_downstream_package/RS_time_up_FDR01_FIXED.txt" "$PKG/01_RNAseq/03_Expression_Tables/GeneList_Increasing_Time_Trend_FDR001.txt"
copy_file "$OUT/rnaseq_downstream_package/RS_time_down_FDR01_FIXED.txt" "$PKG/01_RNAseq/03_Expression_Tables/GeneList_Decreasing_Time_Trend_FDR001.txt"
copy_file "$OUT/rnaseq_downstream_package/RS_all_tested_gene_background_FIXED.txt" "$PKG/01_RNAseq/03_Expression_Tables/Reference_Testable_Genes_Background.txt"

# Search the project for additional ATAC/multiomics artefacts (tables and figures only).
# Do not recopy the already packaged RNA-seq analysis directories or this new package.
while IFS= read -r -d '' f; do
  rel="${f#"$ROOT/"}"
  case "$rel" in
    .git/*|1_data_raw/4_rnaseq/*|reports/"$PACKAGE_NAME"/*) continue ;;
  esac
  lower="${rel,,}"
  [[ "$lower" =~ atac|multiomic|chipseq|multi-omics|differential_peak|peak_annotation ]] || continue
  allowed_ext "$f" || continue
  copy_file "$f" "$PKG/05_ATAC_and_Multiomics/Discovered_Project_Files/$rel"
done < <(find "$ROOT" -type d \( -name .git -o -path "$PKG" \) -prune -o -type f -print0)

# Copy relevant downloaded web-result exports if their filenames identify the tool.
if [[ -d "$HOME/Downloads" ]]; then
  while IFS= read -r -d '' f; do
    allowed_ext "$f" || continue
    base="$(basename "$f")"
    lower="${base,,}"
    if [[ "$lower" =~ panther|gprofiler|g-profiler ]]; then
      copy_file "$f" "$PKG/03_Functional_Enrichment/$(date -r "$f" +%Y%m%d)_$base"
    elif [[ "$lower" =~ lisa|cistrome ]]; then
      copy_file "$f" "$PKG/02_YAP_TEAD_LISA/$(date -r "$f" +%Y%m%d)_$base"
    elif [[ "$lower" =~ wi.?38|telos|rnaseq|rna_seq|alphagenome ]]; then
      copy_file "$f" "$PKG/07_Other_Downloaded_Results/$(date -r "$f" +%Y%m%d)_$base"
    fi
  done < <(find "$HOME/Downloads" -maxdepth 1 -type f -print0)
fi

cat > "$PKG/00_README_FIRST.md" <<'EOF2'
# Telos NGS Project — Presentation Assets

Generated on: @BUILT_AT@

This package collects report-ready figures, result tables, gene lists, selected analysis scripts, project documentation, AlphaGenome outputs, and matching ATAC/multiomics artefacts found in the project. It is organised to make PowerPoint preparation easier.

## Main folders
- `00_Project_Docs`: top-level README and report documents.
- `01_RNAseq`: expression tables, QC, original exploratory plots, downstream tables, and clearly renamed presentation copies of major figures.
- `02_YAP_TEAD_LISA`: prepared increasing/decreasing gene lists and LISA-related outputs found in the project or Downloads.
- `03_Functional_Enrichment`: functional-enrichment exports found in Downloads and other bundled enrichment items.
- `04_AlphaGenome`: AlphaGenome prediction tables, plots, and related output files.
- `05_ATAC_and_Multiomics`: matching ATAC/multiomics files discovered elsewhere in the project.
- `06_Methods_and_Scripts`: analysis scripts and method notes.
- `07_Other_Downloaded_Results`: other relevant downloads found by filename.

## Interpretation cautions
- The workbook's `RS_deseq2` sheet is study-provided; our downstream script did not refit DESeq2.
- Its coefficient represents the supplied time-trend model and is not automatically a direct young-versus-late endpoint test.
- TPM endpoint fold changes are descriptive, not a substitute for count-based differential-expression testing.
- AlphaGenome prediction tracks are model predictions, not observed young/late WI-38 expression measurements.
- LISA nominates candidate regulators using public regulatory datasets; it does not by itself prove binding in our WI-38 samples.
- Recheck custom-background and ID-mapping settings when finalising g:Profiler/PANTHER enrichment.

Raw FASTQ/BAM files, genome indexes, large binary tracks, model weights, and `.git` history are deliberately excluded to keep the presentation ZIP manageable.
EOF2
sed -i "s|@BUILT_AT@|$(date '+%Y-%m-%d %H:%M:%S %Z')|" "$PKG/00_README_FIRST.md"

# If no ATAC results were present, leave an explicit note instead of an empty mystery folder.
if ! find "$PKG/05_ATAC_and_Multiomics/Discovered_Project_Files" -type f -print -quit | grep -q .; then
  cat > "$PKG/05_ATAC_and_Multiomics/README_NO_MATCHING_FILES_FOUND.txt" <<'EOF2'
No project files with ATAC/multiomics-related names and supported result/document formats were discovered automatically in this run.
If teammates have ATAC plots/tables elsewhere, copy them into this folder before presenting.
EOF2
fi

# Build a manifest of package contents and sizes.
find "$PKG" -type f ! -name MANIFEST.tsv -printf '%P\t%s bytes\n' | sort > "$PKG/MANIFEST.tsv"

if ! command -v zip >/dev/null 2>&1; then
  echo "ERROR: 'zip' command is not installed. Install it with: sudo apt install zip" >&2
  echo "The curated folder is still available at: $PKG" >&2
  exit 2
fi
(
  cd "$REPORTS_DIR"
  zip -qr "$ZIP_PATH" "$PACKAGE_NAME"
)

echo
echo "Presentation package folder: $PKG"
echo "ZIP archive:                $ZIP_PATH"
echo "Files included:             $(find "$PKG" -type f | wc -l)"
echo "Package folder size:        $(du -sh "$PKG" | awk '{print $1}')"
echo "ZIP size:                   $(du -sh "$ZIP_PATH" | awk '{print $1}')"
echo
echo "Top-level contents:"
find "$PKG" -mindepth 1 -maxdepth 1 -printf '%f\n' | sort
