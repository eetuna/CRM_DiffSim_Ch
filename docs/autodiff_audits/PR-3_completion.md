# PR-3 Completion Report

## 1) Scope
- Spec: `CRM_v1_0_Differentiable_Equilibrium_Execution_Plan_Frozen.md`
- PR/CP: PR-3 (PyTorch Autograd Backward TIP_POSITION_ONLY)

## 2) Changes Made
- Files modified:
  - `src/CRM_TorchEquilibrium.cpp` (new)
  - `tests/test_crm_tip_position_gradcheck.py` (new)
- New functions/classes:
  - `CRMCatheterModel::EquilibriumConfig`
  - `CRMCatheterModel::crm_equilibrium_forward`
  - `CRMCatheterModel::crm_equilibrium_backward`
  - `CRMCatheterModel::CheckInputs` (internal static helper)
  - `CRMCatheterModel::BuildShootingParams` (internal static helper)
  - `CRMEquilibriumFunction` (Python autograd wrapper)
- Removed/deprecated code: none

## 3) Spec Compliance Checklist
Acceptance criteria:
- PASS: Forward returns `p_tip` with shape `[B, 3]` from PyTorch.
- PASS: Backward returns gradients only for currents `u` and matches `gradcheck` at specified tolerances for 3 operating points.
- PASS: No eliminated Jacobians and no `pseudoInverse()` used in backward.
- PASS: `Li` is accepted as input, enforced `requires_grad=False`, and produces no gradient.

Audit checklist:
- PASS: Autograd backward uses only `A`, `B`, `y_x`, `y_u_direct` with implicit VJP solve (no differentiation through dogleg iterations).
- PASS: `TIP_POSITION_ONLY` path does not compute or save any rotation-related tensors.
- PASS: Input shapes are `[B, N_act, 3]` with runtime check `N_act==1`.
- PASS: `pseudoInverse()` appears nowhere in the new backward path.

## 4) Mathematical / Autodiff Details
- Jacobians used:
  - `A = ∂F_scaled/∂x_scaled` via `CRM_EquilibriumResidualJacobian_A` (MINPACK -> row-major conversion via canonical helper).
  - `B = ∂F_scaled/∂u` via `CRM_EquilibriumResidualJacobian_B` using raw `x_N._u_zc`.
  - `y_x = ∂p_tip/∂x_scaled = (∂p_tip/∂deltau0) * IVALUE_SCALE_DU` using raw `x_N._p_u0`.
  - `y_u_direct = ∂p_tip/∂u` using raw `x_N._p_zc`.
- Linear systems solved:
  - `A^T λ = y_x^T g_p` (3x3) solved with `Eigen::FullPivLU` (direct LU with pivoting).
- MINPACK/solver layout handling: `MinpackJacobianToRowMajor` is used through `CRM_EquilibriumResidualJacobian_A` for `out_fjac` conversion.

## 5) Tests & Evidence
Commands run:
- `python3 -m pytest tests/test_crm_tip_position_gradcheck.py -s`

Results:
- `tests/test_crm_tip_position_gradcheck.py` passed (gradcheck on 3 operating points).

## 6) Determinism & Safety Verification
- Confirmed: no `pseudoInverse()`, no SVD/damping/regularization.
- Confirmed: no differentiation through solver iterations.
- Confirmed: backward recomputation uses only forward-cached tensors (`u`, `Li`, `deltau0`, `ftip`, `localmin`, `residual_norm`, `cfg_equil`).
- Confirmed: non-convergence policy enforced in backward (raises error if `localmin != 0` or `residual_norm > threshold`).

## 7) Go / No-Go Decision
GO to CP-3.
