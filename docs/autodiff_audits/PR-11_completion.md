# PR-11 Completion

## Summary
- Added hybrid analytic VJP contributions for direct eta -> (m_L, n_L) slots in H.
- Preserved FD VJP for all other H dependencies and added a pure-FD reference path.
- Added unit test comparing hybrid vs FD VJP at three stable points.

## Files
- src/CRM_TorchDynamics.cpp
- tests/test_crm_step_vjp_hybrid.py
