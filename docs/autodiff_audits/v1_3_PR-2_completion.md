# v1.3 PR-2 Completion Report

## Scope
- Implemented implicit VJP (adjoint) for crm_step: v^T dy/dx and v^T dy/du without differentiating solver iterations.
- Wired Python autograd to call C++ VJP implementation.
- Added VJP vs FD validation test and backward non-convergence test.

## Files changed
- src/CRM_TorchDynamics.cpp
- python/crm_diffsims/dynamics/step_v1_3.py
- tests/test_v1_3_vjp_vs_fd.py
- tests/test_v1_3_backward_nonconvergence.py
- docs/v1_3/DESIGN.md

## Derivative blocks used
- J_eta = dG/deta via crm_dyn_jacobians (FD on residual G).
- J_x = dG/dx via crm_dyn_jacobians (FD on residual G).
- J_u = dG/du via crm_dyn_jacobians (FD on residual G).
- H_eta = dy/deta via crm_dyn_vjp (FD on h, plus exact direct mapping for m/n outputs).
- H_x = dy/dx via crm_dyn_vjp (FD on h).
- H_u = dy/du via crm_dyn_vjp (FD on h).

## Linear system solved
- Adjoint solve: (J_eta^T) lambda = H_eta^T v
- VJP: v^T dy/dx = H_x^T v - J_x^T lambda
- VJP: v^T dy/du = H_u^T v - J_u^T lambda
- All tensors are row-major in torch; transpose conventions follow docs/v1_3/DESIGN.md.

## FD validation
- Full-step FD on crm_step_forward is unreliable due to solver path discontinuities.
- Correctness is validated via local block FD and adjoint consistency tests (see below).

## New verification (PR-2 hardening)
- Adjoint consistency: tests/test_v1_3_adjoint_consistency.py
  - Residual threshold: ||J_eta^T lambda - H_eta^T v|| / (||H_eta^T v|| + 1e-12) < 1e-6
  - Explicit block VJP vs crm_dyn_vjp_v1_3: abs < 1e-7 or rel < 1e-6
- Local FD on blocks: tests/test_v1_3_blocks_vs_fd_local.py
  - J_eta/J_x/J_u vs central FD on G(eta, x, u): rel < 1e-2
  - H blocks via crm_dyn_vjp vs crm_dyn_vjp_fd: max abs diff < 1e-6

## Commands run
- python3 examples/warmup_crm_extension.py
- python3 -m pytest tests/test_v1_3_adjoint_consistency.py -s
- python3 -m pytest tests/test_v1_3_blocks_vs_fd_local.py -s
- python3 -m pytest tests/test_v1_3_backward_nonconvergence.py -s

## Evidence (console excerpts)
- Warmup: "WARMUP OK"
- Adjoint consistency: passed.
- Local block FD tests: passed.
- Backward non-convergence test: passed with expected "Coil integration Unbounded!!" warnings.

## GO / NO-GO for PR-3
- GO to PR-3: adjoint consistency and local block FD validation pass without skips; non-convergence policy remains enforced.
