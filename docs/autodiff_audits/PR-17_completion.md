# PR-17 Completion

## Summary
- Added outlier drill-down to stress replay test with top-K error indices per dataset and init mode.
- Wrote outlier CSV to output_data/stress_outliers_943.csv with error context, residuals, and classification.
- Added inconsistency classification when large tip error occurs with small residuals.

## Files
- tests/test_stress_replay_npz.py
- docs/autodiff_audits/PR-17_completion.md
- docs/autodiff_audits/CP-17_completion.md
