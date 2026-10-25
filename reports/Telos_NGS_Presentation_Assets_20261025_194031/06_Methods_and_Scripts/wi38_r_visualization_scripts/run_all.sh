#!/usr/bin/env bash
set -euo pipefail

# Run from project root, or pass the project root as the first argument.
if [[ $# -ge 1 ]]; then
  cd "$1"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "$SCRIPT_DIR/01_check_inputs.R"
Rscript "$SCRIPT_DIR/02_volcano_exploratory.R"
Rscript "$SCRIPT_DIR/03_heatmap_top_genes.R"
Rscript "$SCRIPT_DIR/04_pca_samples.R"

echo "All scripts completed. See 1_data_raw/4_rnaseq/outputs/visualizations/"
