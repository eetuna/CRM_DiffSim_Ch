# PR-6 Completion Report

## 1) Scope
- Spec: `Milestone_v1_1_Revised_Differentiable_1_step_dynamics.md`
- PR/CP: PR-6 (eta identification + residual G callable + Jacobian layout contract)

## 2) Changes Made
- Files modified:
  - `src/CRM_TorchDynamics.cpp`
  - `tests/test_crm_dyn_residual.py` (new)
- New functions/classes:
  - `CRMCatheterModel::crm_dyn_residual`
  - `CRMCatheterModel::EvaluateDynResidual` (internal static helper)
- Removed/deprecated code: none

## 3) Spec Compliance Checklist
Acceptance criteria:
- PASS: eta contract implemented in code (`eta = [m_L/IVALUE_SCALE_M, n_L/IVALUE_SCALE_N]` per actuator; dim = `NUM_DYN_RESIDUAL`).
- PASS: G evaluation callable at eta* via `crm_dyn_residual` (no iteration dependence).
- PASS: eta* and residual norm printed/validated in test output.
- PASS: Canonical MINPACK layout conversion helper already present and used for any `fjac` conversions (no ad-hoc indexing).

## 4) Mathematical / Autodiff Details
- Eta definition (code-locked): solver decision vector `eta` is the scaled `[m_L, n_L]` block per actuator, size `6*NUM_ACT_SET`.
- Residual `G` evaluation: uses `DYNNLEquation(eta, ...)` with `DYNNLEqnParams` prepared via `CRMDYNSolverIVP_Prep`, matching the solver call chain.
- No Jacobian matrices solved in PR-6.

## 5) Tests & Evidence
Commands run:
- `python3 -m pytest tests/test_crm_dyn_residual.py -s`

Results (printed by test):
- `eta* tensor([[-3.8328e-01, -7.4667e-04,  2.4701e-04, -3.3566e-05,  2.7008e-03, -7.1174e-01]])`
- `residual_norm tensor([2.5053e-05])`
- Residual norm computed from `G` matches cached `residual_norm` (`torch.allclose` check).

## 6) Determinism & Safety Verification
- Residual evaluation is deterministic and recomputable from cached inputs (`x_t`, `u_t`, `Li`, `eta`, `cfg_dyn`).
- No `pseudoInverse()`, SVD, damping, or regularization introduced in PR-6.

## 7) Go / No-Go Decision
GO to CP-6.
