#!/usr/bin/env bash
set -euo pipefail

STAMP=$(date +%Y%m%d_%H%M%S)
DATASET_PATH=${DATASET_PATH:-$(ls -t output_data/c1_bc_dataset_*.npz | head -1)}

if [[ ! -f "${DATASET_PATH}" ]]; then
  echo "No dataset found at ${DATASET_PATH}" >&2
  exit 1
fi

MODEL_PATH="output_data/models/c1_bc_policy_${STAMP}.pt"
TRAIN_NPZ="output_data/c1_bc_train_${STAMP}.npz"
TRAIN_REPORT="docs/control/run_reports/c1_train_${STAMP}.md"

python3 examples/train_bc_policy.py \
  --dataset "${DATASET_PATH}" \
  --epochs 50 \
  --batch-size 128 \
  --hidden 128,128 \
  --include-tip \
  --smooth-weight 0.1 \
  --checkpoint-every 10 \
  --seed 0 \
  --output-model "${MODEL_PATH}" \
  --output-npz "${TRAIN_NPZ}" \
  --report-path "${TRAIN_REPORT}"

EVAL_PREFIX="c1_bc_eval_${STAMP}"
EVAL_REPORT_PREFIX="c1_eval_${STAMP}"

python3 examples/eval_bc_policy_rollout.py \
  --model "${MODEL_PATH}" \
  --start-idx 0 \
  --subset-len 200 \
  --plot \
  --output-prefix "${EVAL_PREFIX}" \
  --report-prefix "${EVAL_REPORT_PREFIX}" \
  --bc-dataset "${DATASET_PATH}"

