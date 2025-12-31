# CP-14 Completion

## Command
python3 -m pytest tests/test_crm_step_C_step.py -s

## Results
- max_abs=7.248e-05
- max_rel=1.317e-02

## Tolerance
- abs <= 1e-4, rel <= 2e-2
- Justification: C_step_new uses scale-aware steps while FD reference uses fixed step; small deviation is expected from step-size differences while preserving the same residual path.

## Status
GO
