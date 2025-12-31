# PR-13 Completion

## Summary
- Implemented B_step_new with scale-aware u_t steps and param rebuilds in crm_dyn_jacobians.
- Preserved A_step and C_step behavior, and kept FD reference via crm_dyn_jacobians_fd.
- Added unit test comparing B_step_new vs FD reference.

## Files
- src/CRM_TorchDynamics.cpp
- tests/test_crm_step_B_step.py
