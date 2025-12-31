# CP-13 Completion

## Command
python3 -m pytest tests/test_crm_step_B_step.py -s

## Results
- max_abs=0.000e+00
- max_rel=0.000e+00

## Tolerance
- abs <= 1e-5, rel <= 1e-4
- Justification: B_step_new uses scale-aware steps but the same residual evaluation path; close agreement expected.

## Status
GO
