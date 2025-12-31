# PR-15 Completion

## Summary
- Gated FD reference Jacobian/VJP paths behind CRM_DIFFSIMS_ENABLE_FD_REFERENCE.
- Kept v1.2 Jacobians as the default path and preserved hybrid H VJP determinism.
- Added a backward benchmark script comparing FD reference vs v1.2 path.
- Updated tests that rely on FD reference to set the debug flag.

## Files
- src/CRM_TorchDynamics.cpp
- tests/test_crm_step_vjp_hybrid.py
- tests/test_crm_step_A_step_analytic.py
- tests/test_crm_step_B_step.py
- tests/test_crm_step_C_step.py
- examples/benchmark_step_backward.py
- docs/autodiff_audits/PR-15_completion.md
- docs/autodiff_audits/CP-15_completion.md
- docs/autodiff_audits/v1_2_final_completion_summary.md
