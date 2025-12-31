# PR-14 Completion

## Summary
- Implemented C_step_new with scale-aware x_t stepping and full param rebuilds in crm_dyn_jacobians.
- Preserved A_step/B_step behavior and FD reference via crm_dyn_jacobians_fd.
- Added unit test comparing C_step_new vs FD reference.

## Files
- src/CRM_TorchDynamics.cpp
- tests/test_crm_step_C_step.py
