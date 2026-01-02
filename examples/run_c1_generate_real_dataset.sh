#!/usr/bin/env bash
set -euo pipefail

STAMP=$(date +%Y%m%d_%H%M%S)
OUT_NPZ="output_data/c1_bc_dataset_${STAMP}.npz"
REPORT="docs/control/run_reports/c1_generate_dataset_${STAMP}.md"

python3 examples/generate_bc_dataset.py \
  --preset mixed_easy \
  --episodes 400 \
  --episode-len 400 \
  --wrap-dataset \
  --filter-policy drop_unbounded \
  --seed 0 \
  --output-npz "${OUT_NPZ}" \
  --report-path "${REPORT}"

echo "dataset_npz=${OUT_NPZ}"
