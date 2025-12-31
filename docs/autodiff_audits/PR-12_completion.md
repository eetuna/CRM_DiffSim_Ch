# PR-12 Completion

## Summary
- Implemented A_step in crm_dyn_jacobians using scale-aware steps for eta (IVALUE_SCALE_M/N) for NUM_ACT_SET=1.
- Preserved FD B_step/C_step and existing hybrid H VJP.
- Added FD reference API for A_step comparison.

## Files
- src/CRM_TorchDynamics.cpp
- tests/test_crm_step_A_step_analytic.py
