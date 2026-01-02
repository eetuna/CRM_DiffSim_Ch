# Phase C2 MPC-Only Fix

## Root Cause
iLQR could throw an index error when the rollout terminated early, leaving `tip_seq` shorter than `u_seq + 1`, which caused `total_cost` to index past the end.

## Fix
- Truncate `u_seq` to match the rollout length and slice targets accordingly inside `ilqr_solve`.
- Hardened `mpc_only` loop to capture exceptions and report them without killing the overall eval.
- Added diagnostics to the run report for MPC-only failures.

## Commands
```bash
PYTHONPATH=. TORCH_EXTENSIONS_DIR=/workspaces/CRM_DiffSim_Ch/build_torch_ext_dyn \
python3 examples/eval_policy_warmstart_mpc.py \
  --model output_data/models/phase_c2_policy.pt \
  --dataset data/dyn_fk_ramp_circle1_hold1.npz \
  --start-idx 0 \
  --window-len 50 \
  --init-mode dataset_x0 \
  --preroll-steps 5 \
  --horizon 5 \
  --mpc-steps 3 \
  --max-wall-sec 300 \
  --heartbeat-every 1 \
  --seed 0 \
  --output-prefix phase_c2_warmstart_debug \
  --report-path docs/control/run_reports/phase_c2_warmstart_debug.md

pytest -q tests/test_phase_c2_mpc_only_keeps_steps.py -s
```

## Outcome
- `mpc_only` now produces kept steps (`mpc_only_kept: 3`) with finite mean/max errors.
- `mpc_only` exceptions no longer zero out earlier steps; keep masks are preserved.
- Report includes MPC-only exception diagnostics when any iLQR error occurs.

## Evidence (Latest Run)
- Command: `examples/eval_policy_warmstart_mpc.py --model output_data/models/phase_c2_policy.pt --dataset data/dyn_fk_ramp_circle1_hold1.npz --start-idx 0 --window-len 50 --init-mode dataset_x0 --preroll-steps 5 --horizon 5 --mpc-steps 3 --max-wall-sec 300 --heartbeat-every 1 --seed 0 --output-prefix phase_c2_warmstart_debug --report-path docs/control/run_reports/phase_c2_warmstart_debug.md`
- Report: `docs/control/run_reports/phase_c2_warmstart_debug.md`
- Key metrics:
  - `mpc_only_kept: 3`
  - `mpc_only_mean: 10.005735862123784`
  - `mpc_only_max: 15.34742372299618`
  - `eval_go: True`
