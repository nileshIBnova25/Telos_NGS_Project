# Telos NGS Project — Presentation Assets

Generated on: 2026-10-25 19:40:34 IST

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
