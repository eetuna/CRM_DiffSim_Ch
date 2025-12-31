# v1.2 Final Completion Summary

## APIs and Tensor Shapes
- v1.0 equilibrium API (C++ extension):
  - Forward: crm_equilibrium_forward(u, Li, cfg) -> p_tip[B,3], deltau0[B,3], ftip[B,3], localmin[B], residual_norm[B]
  - Backward: crm_equilibrium_backward(grad_p, u, Li, deltau0, ftip, localmin, residual_norm, cfg) -> grad_u[B,1,3]
- v1.1/v1.2 step API (C++ extension):
  - Forward: crm_step_forward(x_t, u_t, Li, cfg) -> x_tp1[B,39], eta[B,6], solver_exit[B], residual_norm[B], u0[B,3], tau[B,1,3]
  - Residual: crm_dyn_residual(eta, x_t, u_t, Li, cfg) -> residual[B,6], u0[B,3], tau[B,1,3]
  - Jacobians: crm_dyn_jacobians(eta, x_t, u_t, Li, cfg) -> A_step[B,6,6], B_step[B,6,3], C_step[B,6,39]
  - VJP: crm_dyn_vjp(grad_x_tp1, eta, x_t, u_t, Li, cfg) -> vjp_eta[B,6], vjp_u[B,1,3], vjp_x[B,39]

## Cached for Backward
- v1.0 equilibrium backward uses cached: u, Li, deltau0, ftip, localmin, residual_norm (from forward), plus cfg.
- v1.1/v1.2 step backward (Python FD bridge) caches: x_t, u_t, Li, eta, solver_exit, residual_norm, cfg.
- v1.1/v1.2 Jacobian/VJP APIs are stateless: they recompute from inputs (eta, x_t, u_t, Li, grad_x_tp1).

## Linear Systems Solved
- v1.0 equilibrium backward solves a 3x3 linear system A^T * lambda = y_x^T * g_p (A = dF/d(deltau0)).
  - Transpose convention: uses A^T for the solve; gradient is y_u^T g_p - B^T lambda.
- v1.1/v1.2 step Jacobian/VJP paths do not solve linear systems; they use residual evaluations and FD or hybrid VJP products.

## Analytic vs FD
- v1.0:
  - Analytic Jacobians for equilibrium residual (CRM_EquilibriumResidualJacobian_A/B) and IVP tip Jacobians (CRMSolverIVPJacobian_TipRaw).
  - Backward uses analytic A/B and raw IVP Jacobians.
- v1.1:
  - Step backward via FD in Python (crm_step_forward perturbations).
  - Jacobians/VJPs provided by FD in crm_dyn_jacobians and crm_dyn_vjp.
- v1.2:
  - A_step/B_step/C_step use scale-aware FD steps (default path).
  - H VJP is hybrid: analytic direct mapping for eta->(m_L,n_L) slots, FD for remaining dependencies.
  - FD reference paths remain available only with CRM_DIFFSIMS_ENABLE_FD_REFERENCE=1.

## Remaining Limitations / v1.3 Ideas
- Full analytic DYNNLEquation Jacobians (A/B/C) are still FD-based; derive analytic residual Jacobians to remove FD cost.
- H VJP still uses FD for nontrivial dependencies; implement analytic adjoints through DYNSolverIVP and CRMFlexible_IVP_Back.
- Add explicit implicit-function backward for crm_step to replace Python FD bridge.
- Expand NUM_ACT_SET > 1 coverage for analytic paths.
