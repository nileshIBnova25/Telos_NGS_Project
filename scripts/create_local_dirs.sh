#!/usr/bin/env bash
set -euo pipefail

# Run from the repository root. These directories contain local data and are ignored by Git.
mkdir -p \
  1_data_raw/0_reference/hg38 \
  1_data_raw/1_preprocessing_and_intermediate \
  1_data_raw/2_benchmark/logs \
  1_data_raw/3_alphagenome/logs \
  1_data_raw/GSE175533_RAW \
  1_data_raw/archive

printf 'Created local TELOS data directories under %s/1_data_raw\n' "$PWD"
printf 'Reminder: files under 1_data_raw/ are intentionally excluded from Git.\n'
