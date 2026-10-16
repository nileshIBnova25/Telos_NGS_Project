# WI-38 Young vs replicative-senescence: AlphaGenome regulatory-locus workflow

This staged workflow reuses the project's existing files and creates a small, auditable AlphaGenome scan around selected differentially expressed genes.

## Input files expected

- `1_data_raw/4_rnaseq/outputs/rnaseq_groundtruth_wi38_young_vs_rs.csv`
- `1_data_raw/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx`, sheet `RS_TC_raw_TPM`
- TPM samples: `PDL20_TP1_A`, `PDL20_TP1_B`, `PDL20_TP1_C`, `PDL46_TP6_A`, `PDL46_TP6_B`, `PDL50_TP7_A`, `PDL50_TP7_B`

The first script uses `python-calamine` because the workbook may be Strict Open XML and may not expose its sheets through `openpyxl`.

## 1. Create a separate Python 3.11 environment

From the repository root:

```bash
conda create -n wi38-alphagenome python=3.11 -y
conda activate wi38-alphagenome
python -m pip install --upgrade pip
python -m pip install pandas numpy python-calamine pyarrow alphagenome plotly
```

AlphaGenome access requires an API key. Request it through the official AlphaGenome site / official Google DeepMind AlphaGenome GitHub documentation. Keep the key private. AlphaGenome API access is intended for non-commercial use under its terms unless separately licensed.

## 2. Put these scripts in the project

Extract this ZIP at the project root. It should create:

`1_data_raw/wi38_r_visualization_scripts/alphagenome_workflow/`

Run commands below from the project root. The scripts use project-relative paths.

```bash
python 1_data_raw/wi38_r_visualization_scripts/alphagenome_workflow/01_prepare_candidates.py
```

Outputs:

- `1_data_raw/4_rnaseq/outputs/alphagenome/candidate_gene_scores_all.csv`
- `1_data_raw/4_rnaseq/outputs/alphagenome/candidate_genes_top20_each_direction.csv`

This chooses up to 20 candidates per direction, with a minimum mean TPM of 1 in at least one group. It prioritizes agreement between ground-truth log2FC and the TPM group means, plus replicate-pattern consistency. Inspect the CSV before proceeding; do not blindly interpret a large fold change as robust if TPM is low or replicates disagree.

## 3. Set API key and inspect track metadata

In the same terminal session:

```bash
read -rsp 'AlphaGenome API key: ' ALPHAGENOME_API_KEY; echo
export ALPHAGENOME_API_KEY
python 1_data_raw/wi38_r_visualization_scripts/alphagenome_workflow/02_inspect_alphagenome_tracks.py
```

The script writes full AlphaGenome human-output metadata and a fibroblast-related subset into `1_data_raw/4_rnaseq/outputs/alphagenome/`. Inspect `alphagenome_candidate_ontology_terms.csv` and the metadata CSV. Choose ontology CURIEs for relevant fibroblast or closely related biosamples only if they actually appear. Do not assume any track is WI-38-specific or represents replicative senescence. If no appropriate fibroblast track exists, document the closest available reference biosample as a limitation.

## 4. Run a small pilot regulatory scan

After reviewing the actual ontology CURIEs, provide one or more exact IDs from the metadata. Do not use a term merely because it sounds relevant.

```bash
export ALPHAGENOME_ONTOLOGY_TERMS='REPLACE_WITH_REVIEWED_CURIE_1,REPLACE_WITH_REVIEWED_CURIE_2'
# Start with five up and five down candidates (default); optionally change count:
export ALPHAGENOME_N_GENES_PER_DIRECTION=5
python 1_data_raw/wi38_r_visualization_scripts/alphagenome_workflow/03_predict_regulatory_loci.py
```

The prediction script uses GENCODE v46 hg38 annotations and the official AlphaGenome API to predict RNA-seq, CAGE, DNase, ATAC and histone-mark tracks in a 1 Mb reference-genome interval around selected genes. It saves per-gene interactive HTML tracks to `1_data_raw/4_rnaseq/outputs/alphagenome/locus_plots/` and regional signal summaries to `alphagenome_regulatory_signal_summaries.csv`.

## How to interpret

- `RNA_SEQ`: predicted RNA-seq-like signal, not the same measurement as this experiment's TPM.
- `CAGE`: predicted transcription-start-site activity.
- `DNASE` / `ATAC`: predicted chromatin accessibility.
- `CHIP_HISTONE`: predicted histone-mark signal; interpretation depends on the actual track/assay name.
- Compare promoter and near-TSS signal across selected loci and among available biosample tracks. Do not directly subtract arbitrary track values across assays because scales and normalization differ.
- AlphaGenome's reference sequence is generally the same hg38 reference for every WI-38 sample, so these predictions do not independently model the young-to-senescent transition. Regulatory predictions suggest hypotheses to test; they do not establish the causal mechanism of senescence-associated expression changes.
- The candidate-list step is only shortlist/QC. Use a count-based differential-expression method with appropriate replicates and design for formal DE claims where raw counts are available.

## Troubleshooting

- If candidate selection reports low match rate, compare gene IDs in the CSV and TPM workbook columns. Do not hand-map gene IDs to symbols without verifying the mapping.
- If AlphaGenome returns no tracks for a gene, verify that the chosen ontology terms exist for the requested output types in the full metadata. An empty result is not proof that a locus has no activity.
- If `get_gene_interval` cannot map a gene, verify it is an HGNC symbol or Ensembl gene ID compatible with the hg38 GENCODE annotation.
