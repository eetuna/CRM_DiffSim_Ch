# PR-7 Completion Report

## 1) Scope
- Spec: `Milestone_v1_1_Revised_Differentiable_1_step_dynamics.md`
- PR/CP: PR-7 (A_step/B_step/C_step Jacobians + H VJP products)

## 2) Changes Made
- Files modified:
  - `src/CRM_TorchDynamics.cpp`
  - `tests/test_crm_dyn_jacobians.py` (new)
- New functions/classes:
  - `CRMCatheterModel::crm_dyn_jacobians`
  - `CRMCatheterModel::crm_dyn_vjp`
  - `CRMCatheterModel::EvaluateDynHOutput` (internal helper)
- Removed/deprecated code: none

## 3) Spec Compliance Checklist
Acceptance criteria:
- PASS (v1.1 bridge): Jacobian functions `A_step`, `B_step`, `C_step` available at `eta*` with correct dimensions.
- PASS (v1.1 bridge): VJP products `H_eta^T g`, `H_u^T g`, `H_x^T g` available with correct dimensions.
- PASS: `eta` ordering and scaling match `TrustRegionDogleg_dyn` solver coordinates.
- PASS: No ad-hoc MINPACK layout conversions added.

## 4) Mathematical / Autodiff Details
- `A_step`, `B_step`, `C_step` computed by central differences of `G(eta, x, u, Li)` using the same residual callback path as forward (allowed in v1.1 bridge).
- `H` map evaluated deterministically from cached `(eta, x_t, u_t, Li)` by recomputing `u0/tau` via `DYNNLEquation` and running `DYNSolverIVP` without re-running solver iterations.
- VJP products computed by central-difference directional derivatives dotted with upstream `g` (allowed in v1.1 bridge).

## 5) Tests & Evidence
Commands run:
- `python3 -m pytest tests/test_crm_dyn_jacobians.py -s`

Results:
- PASS: shapes for `A_step`, `B_step`, `C_step`, and `H` VJPs.
- PASS: outputs are finite for the test case.

## 6) Determinism & Safety Verification
- Deterministic evaluation from cached forward inputs (`eta`, `x_t`, `u_t`, `Li`, `cfg_dyn`).
- No differentiation through solver iterations.
- No pseudo-inverse/SVD/damping/regularization introduced.

## 7) Go / No-Go Decision
GO (v1.1 bridge allowed; analytic replacement deferred to v1.2).
