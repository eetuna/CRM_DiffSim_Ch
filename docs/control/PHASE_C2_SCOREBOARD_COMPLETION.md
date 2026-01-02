# Phase C2 Scoreboard Completion

## Summary

- Added `examples/scoreboard_eval.py` and smoke tests for scoreboard + CEM dataset generation.
- Switched dataset generation to accept `--controller cem_mpc` with safe-bounds overrides and explicit keep/drop masks.
- Generated scoreboard artifacts for replay, policy, and policy+MPC on the circle window.

## Command Outputs

### 1) Warmup extension

```
TORCH_EXTENSIONS_DIR=/workspaces/CRM_DiffSim_Ch/build_torch_ext_dyn
compiling/ext-loading... elapsed_sec=0.0
WARMUP OK
```

### 2) Scoreboard smoke test

```
.
1 passed in 6.24s
```

### 3) CEM dataset smoke test

```
Saved dataset to output_data/test_bc_dataset_cem_smoke.npz
dataset_counts: total=20 kept=20 dropped=0 reasons={"unbounded": 0, "solver_exit": 0, "residual": 0, "nonfinite": 0, "target_nonfinite": 0, "policy_drop": 0, "prefix_fail": 0}
Run report: docs/control/run_reports/test_bc_dataset_cem_smoke.md
.
1 passed in 7.77s
```

### 4) Scoreboard eval (50-step window)

```
TORCH_EXTENSIONS_DIR=/workspaces/CRM_DiffSim_Ch/build_torch_ext_dyn
step_start mode=replay_only step=0
step_end mode=replay_only step=0 step_sec=0.19
step_start mode=replay_only step=1
step_end mode=replay_only step=1 step_sec=0.01
step_start mode=replay_only step=2
step_end mode=replay_only step=2 step_sec=0.00
step_start mode=replay_only step=3
step_end mode=replay_only step=3 step_sec=0.01
step_start mode=replay_only step=4
step_end mode=replay_only step=4 step_sec=0.01
step_start mode=replay_only step=5
step_end mode=replay_only step=5 step_sec=0.01
step_start mode=replay_only step=6
step_end mode=replay_only step=6 step_sec=0.01
step_start mode=replay_only step=7
step_end mode=replay_only step=7 step_sec=0.01
step_start mode=replay_only step=8
step_end mode=replay_only step=8 step_sec=0.01
step_start mode=replay_only step=9
step_end mode=replay_only step=9 step_sec=0.01
step_start mode=replay_only step=10
step_end mode=replay_only step=10 step_sec=0.01
step_start mode=replay_only step=11
step_end mode=replay_only step=11 step_sec=0.01
step_start mode=replay_only step=12
step_end mode=replay_only step=12 step_sec=0.01
step_start mode=replay_only step=13
step_end mode=replay_only step=13 step_sec=0.01
step_start mode=replay_only step=14
step_end mode=replay_only step=14 step_sec=0.01
step_start mode=replay_only step=15
step_end mode=replay_only step=15 step_sec=0.01
step_start mode=replay_only step=16
step_end mode=replay_only step=16 step_sec=0.00
step_start mode=replay_only step=17
step_end mode=replay_only step=17 step_sec=0.00
step_start mode=replay_only step=18
step_end mode=replay_only step=18 step_sec=0.01
step_start mode=replay_only step=19
step_end mode=replay_only step=19 step_sec=0.01
step_start mode=replay_only step=20
step_end mode=replay_only step=20 step_sec=0.00
step_start mode=replay_only step=21
step_end mode=replay_only step=21 step_sec=0.00
step_start mode=replay_only step=22
step_end mode=replay_only step=22 step_sec=0.01
step_start mode=replay_only step=23
step_end mode=replay_only step=23 step_sec=0.00
step_start mode=replay_only step=24
step_end mode=replay_only step=24 step_sec=0.00
step_start mode=replay_only step=25
step_end mode=replay_only step=25 step_sec=0.00
step_start mode=replay_only step=26
step_end mode=replay_only step=26 step_sec=0.01
step_start mode=replay_only step=27
step_end mode=replay_only step=27 step_sec=0.00
step_start mode=replay_only step=28
step_end mode=replay_only step=28 step_sec=0.00
step_start mode=replay_only step=29
step_end mode=replay_only step=29 step_sec=0.00
step_start mode=replay_only step=30
step_end mode=replay_only step=30 step_sec=0.00
step_start mode=replay_only step=31
step_end mode=replay_only step=31 step_sec=0.01
step_start mode=replay_only step=32
step_end mode=replay_only step=32 step_sec=0.00
step_start mode=replay_only step=33
step_end mode=replay_only step=33 step_sec=0.00
step_start mode=replay_only step=34
step_end mode=replay_only step=34 step_sec=0.00
step_start mode=replay_only step=35
step_end mode=replay_only step=35 step_sec=0.00
step_start mode=replay_only step=36
step_end mode=replay_only step=36 step_sec=0.01
step_start mode=replay_only step=37
step_end mode=replay_only step=37 step_sec=0.00
step_start mode=replay_only step=38
step_end mode=replay_only step=38 step_sec=0.00
step_start mode=replay_only step=39
step_end mode=replay_only step=39 step_sec=0.00
step_start mode=replay_only step=40
step_end mode=replay_only step=40 step_sec=0.01
step_start mode=replay_only step=41
step_end mode=replay_only step=41 step_sec=0.00
step_start mode=replay_only step=42
step_end mode=replay_only step=42 step_sec=0.00
step_start mode=replay_only step=43
step_end mode=replay_only step=43 step_sec=0.00
step_start mode=replay_only step=44
step_end mode=replay_only step=44 step_sec=0.01
step_start mode=replay_only step=45
step_end mode=replay_only step=45 step_sec=0.00
step_start mode=replay_only step=46
step_end mode=replay_only step=46 step_sec=0.00
step_start mode=replay_only step=47
step_end mode=replay_only step=47 step_sec=0.01
step_start mode=replay_only step=48
step_end mode=replay_only step=48 step_sec=0.01
step_start mode=replay_only step=49
step_end mode=replay_only step=49 step_sec=0.01
summary mode=replay_only kept=50/50 mean_mm=5.134 max_mm=39.296 unbounded=0 solver_exit_counts={"0": 50} runtime_sec=0.50
artifact_npz=output_data/scoreboard_20260102_225133_replay_only.npz
report_path=docs/control/run_reports/scoreboard_20260102_225133_replay_only.md
step_start mode=policy_only step=0
step_end mode=policy_only step=0 step_sec=0.02
step_start mode=policy_only step=1
step_end mode=policy_only step=1 step_sec=0.02
step_start mode=policy_only step=2
step_end mode=policy_only step=2 step_sec=0.01
step_start mode=policy_only step=3
step_end mode=policy_only step=3 step_sec=0.01
step_start mode=policy_only step=4
step_end mode=policy_only step=4 step_sec=0.01
step_start mode=policy_only step=5
step_end mode=policy_only step=5 step_sec=0.01
step_start mode=policy_only step=6
step_end mode=policy_only step=6 step_sec=0.00
step_start mode=policy_only step=7
step_end mode=policy_only step=7 step_sec=0.01
step_start mode=policy_only step=8
step_end mode=policy_only step=8 step_sec=0.01
step_start mode=policy_only step=9
step_end mode=policy_only step=9 step_sec=0.01
step_start mode=policy_only step=10
step_end mode=policy_only step=10 step_sec=0.01
step_start mode=policy_only step=11
step_end mode=policy_only step=11 step_sec=0.01
step_start mode=policy_only step=12
step_end mode=policy_only step=12 step_sec=0.01
step_start mode=policy_only step=13
step_end mode=policy_only step=13 step_sec=0.02
step_start mode=policy_only step=14
step_end mode=policy_only step=14 step_sec=0.02
summary mode=policy_only kept=14/15 mean_mm=16.736 max_mm=39.296 unbounded=0 solver_exit_counts={"0": 14, "3": 1} runtime_sec=0.17
artifact_npz=output_data/scoreboard_20260102_225133_policy_only.npz
report_path=docs/control/run_reports/scoreboard_20260102_225133_policy_only.md
/workspaces/CRM_DiffSim_Ch/examples/scoreboard_eval.py:193: UserWarning: To copy construct from a tensor, it is recommended to use sourceTensor.detach().clone() or sourceTensor.detach().clone().requires_grad_(True), rather than torch.tensor(sourceTensor).
  inputs = [state, torch.tensor(target_t, dtype=torch.float32)]
iLQR iter 1: cost=5.815085e+03 reg=1.000e-03 clamp_rate=0.00% sign_flip_rate=0.00%
step_start mode=policy_plus_mpc step=0
step_end mode=policy_plus_mpc step=0 step_sec=0.01
/workspaces/CRM_DiffSim_Ch/examples/scoreboard_eval.py:193: UserWarning: To copy construct from a tensor, it is recommended to use sourceTensor.detach().clone() or sourceTensor.detach().clone().requires_grad_(True), rather than torch.tensor(sourceTensor).
  inputs = [state, torch.tensor(target_t, dtype=torch.float32)]
iLQR iter 1: cost=1.437142e+03 reg=5.000e-05 clamp_rate=0.00% sign_flip_rate=0.00%
step_start mode=policy_plus_mpc step=1
step_end mode=policy_plus_mpc step=1 step_sec=0.01
iLQR iter 1: cost=2.523554e+03 reg=1.000e-03 clamp_rate=0.00% sign_flip_rate=0.00%
summary mode=policy_plus_mpc kept=2/2 mean_mm=38.233 max_mm=39.296 unbounded=0 solver_exit_counts={"0": 2} runtime_sec=69.44
artifact_npz=output_data/scoreboard_20260102_225133_policy_plus_mpc.npz
report_path=docs/control/run_reports/scoreboard_20260102_225133_policy_plus_mpc.md
```

## Artifacts

- `output_data/scoreboard_20260102_225133_replay_only.npz`
- `output_data/scoreboard_20260102_225133_policy_only.npz`
- `output_data/scoreboard_20260102_225133_policy_plus_mpc.npz`
- `docs/control/run_reports/scoreboard_20260102_225133_replay_only.md`
- `docs/control/run_reports/scoreboard_20260102_225133_policy_only.md`
- `docs/control/run_reports/scoreboard_20260102_225133_policy_plus_mpc.md`
- `output_data/test_bc_dataset_cem_smoke.npz`
- `docs/control/run_reports/test_bc_dataset_cem_smoke.md`
