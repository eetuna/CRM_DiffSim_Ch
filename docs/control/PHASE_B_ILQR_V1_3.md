# Phase B: iLQR/MPC with v1.3 Linearization

## Overview
This phase runs iLQR/MPC using the v1.3 implicit VJP linearization backend. The controller uses the same cost and state conventions as Phase A and enforces safe regime clamps (|u|, |Δu|) from `docs/control/safe_bounds.json` when available.

## How to run
Circle:

```
python3 examples/run_ilqr_circle.py --backend v1_3 --plot
```

Lemniscate:

```
python3 examples/run_ilqr_lemniscate.py --backend v1_3 --plot
```

Lemniscate (staged tracking recommended):

```
python3 examples/run_ilqr_lemniscate.py --backend v1_3 --mode staged --plot
```

Notes for v1.3 example defaults:
- The v1.3 examples use reduced horizons/iterations to keep runtime bounded.

## Failure modes and guardrails
- If forward nonconvergence is detected (solver exit != 0 or residual > threshold), the rollout is treated as invalid. The line search backs off and regularization increases.
- Clamp enforcement:
  - |u| <= umax_safe
  - |Δu| <= d_umax_safe
  - sign-flip caps apply to prevent instantaneous flips.
- `reg` (λ) controls trust region strength; larger values penalize aggressive updates.
- Line search alpha is reduced on failed rollouts or cost regressions.

## Lemniscate modes: hold vs staged vs track
- `hold`: hold the initial tip target for the entire run; safest but does not track the lemniscate.
- `staged`: hold for `N_warm` steps, then blend into the lemniscate over `N_blend` steps; recommended default for full trajectory tracking.
- `track`: use the lemniscate targets immediately; can be more sensitive to rollout instability.

## Plots and inspection
Use `--plot` to save PNG+PDF plots to `output_data/` (tip trajectory, currents, error vs time). If matplotlib is unavailable, plotting is skipped with a warning.

## Notes
- v1.3 backend uses operator linearization internally and never calls autograd Jacobian on the dynamics step.
- Keep all runs inside the safe operating regime for stability.
