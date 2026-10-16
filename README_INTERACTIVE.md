# Interactive WI-38 RNA-seq plots (HTML)

This adds hover-enabled HTML versions of the volcano plot, heatmap, and PCA. Hover a volcano point to see the gene and statistics; hover a heatmap cell to see the gene, sample, TPM, ground-truth log2FC, and row Z-score; hover a PCA point to see the sample and PC scores.

## Run from your project root

The project root is the folder containing `1_data_raw/`:

```bash
cd ~/Sem3/BDBP307_Omics_Lab/ProjectNGS/Telos_NGS_Project
python3 -m pip install --user pandas numpy openpyxl plotly
python3 1_data_raw/wi38_r_visualization_scripts/05_interactive_html.py
```

If you use a virtual environment, activate it first and omit `--user` from the pip command.

The script expects these files:

- `1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx` (sheet `RS_TC_raw_TPM`)
- `1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv`
- `1_data_raw/4_rnaseq/outputs/visualizations/exploratory_limma_results.csv` (created by the existing R volcano script)

It writes to `1_data_raw/4_rnaseq/outputs/visualizations/interactive/`. Open `index.html` in Firefox or Chrome. Keep `plotly.min.js` in the same folder as the HTML files for offline viewing.

## Interpreting the figures carefully

- **PCA:** The supplied static plot has the young PDL20 replicates clustered closely, whereas the late/RS samples are displaced from them. The two PDL46/TP6 samples also separate from the two PDL50/TP7 samples. This could reflect biological progression, passage/timepoint effects, sample variation, or other confounding; with only seven samples, these possibilities cannot be disentangled from PCA alone.
- **Heatmap:** Colors are row-wise Z-scores of `log2(TPM + 0.001)`. A red cell means higher expression relative to that gene's other samples, not necessarily high absolute TPM. Hover cells to inspect the measurements.
- **Volcano:** Positive log2FC means higher in RS; negative means lower in RS. The displayed p-values came from exploratory limma on transformed TPM, not a count-based differential-expression model. Do not report these p-values as definitive RNA-seq differential-expression evidence. Formal analysis should use appropriate counts/estimated counts and a method such as DESeq2, edgeR, or limma-voom with a design matching the experiment.

Any candidate genes suggested by the volcano plot should be checked against their per-sample TPM values, gene identifiers, and biological context before being interpreted as senescence or myofibroblast markers.
