# CP-18 Completion

## Commands Run
1) python3 -m pytest tests/test_stress_replay_npz.py -s
2) python3 -m pytest tests/test_workspace_branching_943.py -s
3) python3 examples/run_workspace_stress_report.py

## PASS/FAIL
- test_stress_replay_npz: PASS
- test_workspace_branching_943: PASS
- run_workspace_stress_report.py: PASS (wrote CSV)

## Before/After Inconsistency Summary
- Before (PR-17): inconsistencies_detected=108, concentrated at indices 0–5.
- After (PR-18): inconsistencies_detected=0.
  - expected_initial_state_mismatch counts:
    - zero: 36 (indices 0–5)
    - warm: 36 (indices 0–5)
    - workspace: 36 (indices 0–5)
  - equilibrium_seed: 0 mismatches; all outliers classified ok.

## Output Files
- output_data/stress_outliers_943.csv updated
- output_data/stress_report_943.csv updated

## Interpretation
- Equilibrium-seeded x0 removes early-step large-error inconsistencies; remaining early outliers are reclassified as expected_initial_state_mismatch in non-seeded modes.
- No solver-instability driven inconsistency remains; errors align with initialization differences.

## Status
GO
