# Safe Operating Regime (Forward Sweep)

Safe bounds are derived empirically from the captured sweep CSV:
- `output_data/forward_stability_sweep_captured.csv`

## Classification

- SAFE: `unbounded_detected=False` and `nan_detected=False`
- WARNING: `unbounded_detected=True` and `nan_detected=False`
- FAIL: `nan_detected=True`

## Summary table

| umax | d_umax | status |
| --- | --- | --- |
| 0.05 | 0.01 | SAFE |
| 0.05 | 0.02 | SAFE |
| 0.10 | 0.01 | SAFE |
| 0.10 | 0.02 | WARNING |
| 0.20 | 0.01 | WARNING |
| 0.20 | 0.02 | SAFE |
| 0.30 | 0.01 | WARNING |
| 0.30 | 0.02 | WARNING |

## Recommended default

- umax_safe = 0.10
- d_umax_safe = 0.01

These bounds are empirical, derived from forward-only stability sweeps, and should be enforced by all controllers/optimizers until new sweep data justifies updates.
