#!/usr/bin/env bash
set -euo pipefail

# TELOS directory setup.
# Run from the repository root.
# This script creates directories only; it does not download,
# move, delete, or process project data.

# ------------------------------------------------------------
# 1. Local data, reference genome and benchmark
# ------------------------------------------------------------
mkdir -p \
  1_data_raw/0_reference/hg38 \
  1_data_raw/1_preprocessing_and_intermediate \
  1_data_raw/2_benchmark/logs \
  1_data_raw/3_alphagenome/00_pilot/logs \
  1_data_raw/3_alphagenome/01_config \
  1_data_raw/3_alphagenome/02_atac_run/logs \
  1_data_raw/4_rnaseq/inputs \
  1_data_raw/4_rnaseq/outputs \
  1_data_raw/4_rnaseq/qc \
  1_data_raw/4_rnaseq/logs \
  1_data_raw/GSE175533_RAW \
  1_data_raw/archive

# ------------------------------------------------------------
# 2. Analysis scripts
# ------------------------------------------------------------
mkdir -p \
  2_scripts/1_Preprocessing \
  2_scripts/2_AlphaGenome \
  2_scripts/3_Evaluation \
  2_scripts/4_RNAseq \
  2_scripts/5_Integration

# ------------------------------------------------------------
# 3. Derived results
# ------------------------------------------------------------
mkdir -p \
  3_results/tables \
  3_results/figures \
  3_results/reports

# ------------------------------------------------------------
# 4. Inform the user
# ------------------------------------------------------------
printf 'TELOS directory setup validated for repository: %s\n' "$PWD"
printf 'Local data directory: %s/1_data_raw\n' "$PWD"
printf 'Derived results directory: %s/3_results\n' "$PWD"
printf 'Raw data and generated results must be excluded from Git.\n'
printf 'This script creates directories only; it does not process data.\n'
