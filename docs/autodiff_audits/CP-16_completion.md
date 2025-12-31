# CP-16 Completion

## Commands
1) python3 -m pytest tests/test_stress_replay_npz.py -s
2) python3 -m pytest tests/test_workspace_branching_943.py -s

## Outputs Summary
- test_stress_replay_npz: PASS (runtime ~2m17s)
  - conv_rate ~0.994-0.997 across datasets
  - max_residual_norm ~2.25e-02 to 3.97e-02
  - mean_tip_error ~2.64e-01 to 5.63e-01
  - max_tip_error ~3.44e+01
  - best_init reported as "zero" (ties across modes)
- test_workspace_branching_943: PASS
  - sign_agreement_rate=0.9766
  - positive_y_rate=0.9996
  - negative_y_rate=0.9547
  - ambiguous_fraction=0.0000

## Interpretation
- Stability: high convergence rate across all datasets with residuals below the spike threshold; no catastrophic divergence.
- Branching: strong sign agreement between tip_y and c3 with low ambiguity in near-zero region.
- Initialization behavior: zero, warm, and workspace initializations show similar convergence and error profiles on these datasets.

## Status
GO
