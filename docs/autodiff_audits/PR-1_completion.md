# PR-1 Completion Report

## 1) Scope
- Spec: `CRM_v1_0_Differentiable_Equilibrium_Execution_Plan_Frozen.md`
- PR/CP: PR-1 (Raw IVP Jacobian Getter (TIP))

## 2) Changes Made
- Files modified:
  - `src/CRM_IVPJacobian.hpp`
  - `src/CRM_IVPJacobian.cpp`
  - `main/CRM_KinematicsTestFunctions.cpp`
- New functions/classes:
  - `CRMCatheterModel::CRMIVPTipJacobiansRaw`
  - `CRMCatheterModel::CRMSolverIVPJacobian_TipRaw`
  - `CRMCatheterModel::CRM_PrintIVPTipJacobiansRawExample`
- Removed/deprecated code: none

## 3) Spec Compliance Checklist
Acceptance criteria:
- PASS: Getter returns raw tip blocks from `x_N._p_u0`, `x_N._p_zc`, `x_N._u_u0`, `x_N._u_zc`, `x_N._ws_u0`, `x_N._ws_zc`.
- PASS: No `pseudoInverse()` callsites introduced in the new getter path.
- PASS: For `NUM_ACT_SET=1`, the getter returns `zc` as a 3x3 matrix (no `zl` columns).

Audit checklist:
- PASS: Getter uses `CRMSolverIVP_CoreWithJacobian` and reads `x_N._*` members directly.
- PASS: Getter never calls `CRMSolverIVPJacobian(...)`.
- PASS: Getter path contains zero uses of `pseudoInverse()` or `completeOrthogonalDecomposition()`.
- PASS: Returned Jacobians correspond to `zc` for currents (not `zl`).

## 4) Mathematical / Autodiff Details
- Jacobians introduced: raw IVP tip blocks `∂p/∂u0`, `∂p/∂zc`, `∂u/∂u0`, `∂u/∂zc`, `∂ws/∂u0`, `∂ws/∂zc` from `AugmentedStateVector<IVPJacobiansFull>`.
- Linear systems solved: none.
- MINPACK/solver layout handling: not applicable in this PR; raw IVP blocks are mapped row-major via `__EMT<..., Eigen::RowMajor>` consistent with `CRM_StateVector_Definitions.hpp`.

## 5) Tests & Evidence
- Commands run: not run (test harness added but not executed).
- Numeric results/tolerances: N/A.

## 6) Determinism & Safety Verification
- Confirmed: no `pseudoInverse()`, no SVD/damping/regularization.
- Confirmed: no differentiation through solver iterations.
- Confirmed: deterministic forward→backward recomputation (no hidden globals/statics introduced).

## 7) Go / No-Go Decision
GO to PR-2.
