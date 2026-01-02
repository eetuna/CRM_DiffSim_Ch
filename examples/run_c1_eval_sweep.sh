#!/usr/bin/env bash
set -euo pipefail

MODEL_PATH=${MODEL_PATH:-}
if [[ -z "${MODEL_PATH}" ]]; then
  MODEL_PATH=$(ls -t output_data/models/c1_bc_policy_*.pt 2>/dev/null | head -1 || true)
fi
if [[ -z "${MODEL_PATH}" ]]; then
  echo "Set MODEL_PATH to the trained BC model (e.g. output_data/models/c1_bc_policy_*.pt)" >&2
  exit 1
fi
echo "Using MODEL_PATH=${MODEL_PATH}"

STAMP=$(date +%Y%m%d_%H%M%S)
STARTS=(0 100)
WINDOW=80

OUT_PREFIX="c1_eval_sweep_${STAMP}"

for start in "${STARTS[@]}"; do
  for preroll in 0 10; do
    python3 examples/eval_bc_policy_rollout.py \
      --model "${MODEL_PATH}" \
      --init-mode dataset_x0 \
      --circle-dataset data/dyn_fk_ramp_circle1_hold1.npz \
      --lemniscate-dataset data/dyn_fk_lem1_y40_a10_L94_hold1.npz \
      --start-idx "${start}" \
      --subset-len "${WINDOW}" \
      --preroll-steps "${preroll}" \
      --plot \
      --output-prefix "${OUT_PREFIX}_preroll${preroll}_${start}" \
      --report-prefix "c1_eval_sweep_preroll${preroll}_${start}"
  done
done

python3 examples/summarize_bc_eval_npz.py "output_data/${OUT_PREFIX}_*.npz" --plot

python3 - <<'PY'
import glob
import os
import numpy as np

prefix = "output_data/c1_eval_sweep_*_preroll"
paths = glob.glob(f"{prefix}0_150_lemniscate.npz") + glob.glob(f"{prefix}10_150_lemniscate.npz")
payloads = {}
for path in paths:
    data = np.load(path, allow_pickle=True)
    preroll = int(data.get("preroll_steps", 0))
    kept = int(data.get("kept_steps_policy", data.get("policy_keep_mask", np.zeros((0,))).sum()))
    payloads[preroll] = (os.path.basename(path), kept)
for preroll in [0, 10]:
    if preroll in payloads:
        name, kept = payloads[preroll]
        print(f"preroll_steps={preroll} file={name} kept_steps={kept}")
if 0 in payloads and 10 in payloads:
    if payloads[0][1] == 0 and payloads[10][1] > 0:
        print("preroll_steps=10 fixes kept_steps=0 at start_idx=150")
    else:
        print("preroll_steps=10 did not fix kept_steps=0 at start_idx=150")
PY
