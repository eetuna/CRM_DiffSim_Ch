# CP-17 Completion

## Command
python3 -m pytest tests/test_stress_replay_npz.py -s

## Outputs Summary
- test_stress_replay_npz: PASS
- outlier CSV written to output_data/stress_outliers_943.csv
- inconsistencies_detected=108 (large tip error with small residual)

## Interpretation
- Large tip-error outliers are present despite low residuals; flagged as inconsistencies for follow-up.
- Nonconvergence cases are classified as expected where residuals are high.

## Status
GO
