# CP-12 Completion

## Command
python3 -m pytest tests/test_crm_step_A_step_analytic.py -s

## Results
- max_abs=0.000e+00
- max_rel=0.000e+00

## Tolerance
- abs <= 1e-6, rel <= 1e-5
- Justification: A_step uses the same residual evaluations with scale-aware steps; only floating-point noise expected.

## Status
GO
