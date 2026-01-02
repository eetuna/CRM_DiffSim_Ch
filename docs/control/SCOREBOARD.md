# Scoreboard Evaluation

This document describes the lightweight scoreboard evaluation used to compare rollout
quality across multiple controllers on a fixed dataset window.

## Script

Run:

```
PYTHONPATH=. TORCH_EXTENSIONS_DIR=/workspaces/CRM_DiffSim_Ch/build_torch_ext_dyn \
python3 examples/scoreboard_eval.py \
  --dataset data/dyn_fk_ramp_circle1_hold1.npz \
  --start-idx 0 \
  --window-len 50 \
  --modes replay_only,policy_only,policy_plus_mpc,cem_mpc \
  --max-wall-sec 60
```

Outputs are written per mode:
- `output_data/scoreboard_<mode>_<timestamp>.npz`
- `docs/control/run_reports/scoreboard_<mode>_<timestamp>.md`

## Modes

- `replay_only`: open-loop replay of dataset currents with safe bounds clamping.
- `policy_only`: BC policy closed-loop tracking with bounded currents.
- `policy_plus_mpc`: iLQR MPC seeded by the policy warm-start.
- `cem_mpc`: Jacobian-free CEM MPC (open-loop evaluation of the planned sequence).

All modes share the same initialization: dataset `x0` with optional preroll (dataset
controls applied open-loop before the window).

## Metrics

The console summary for each mode includes:
- `kept/attempted`: number of steps passing solver/finite checks.
- `mean_mm`, `max_mm`: mean and max tip error over kept steps (mm).
- `unbounded`: count of steps flagged nonfinite.
- `solver_exit_counts`: solver exit code histogram.
- `runtime_sec`: wall-clock runtime per mode.

Each NPZ stores numeric-only arrays:
- `tip_xyz`, `target_xyz`, `u_t`, `u_raw`, `keep_mask`, `error_mm`
- solver telemetry (`solver_exit`, `residual_norm`, `unbounded_detected`)
- metadata scalars (`dt`, `Li`, `umax`, `d_umax`)

Run reports record the same metrics plus the invocation and dataset slice.
