# PR-2 Completion Report

## 1) Scope
- Spec: `CRM_v1_0_Differentiable_Equilibrium_Execution_Plan_Frozen.md`
- PR/CP: PR-2 (Residual Jacobians A/B + FD test)

## 2) Changes Made
- Files modified:
  - `src/CRM_BVPIVP_APIDeclarations.hpp`
  - `src/CRM_BVPSolver.cpp`
  - `main/CRM_KinematicsTestFunctions.cpp`
- New functions/classes:
  - `CRMCatheterModel::MinpackJacobianToRowMajor`
  - `CRMCatheterModel::CRM_EquilibriumResidual_FScaled`
  - `CRMCatheterModel::CRM_EquilibriumResidualJacobian_A`
  - `CRMCatheterModel::CRM_EquilibriumResidualJacobian_B`
  - `CRMCatheterModel::CRM_FDTest_EquilibriumResidualJacobian_B`
  - `CRMCatheterModel::CRM_PrepareNLEParamsFromShooting` (internal static helper in `src/CRM_BVPSolver.cpp`)
- Removed/deprecated code: none

## 3) Spec Compliance Checklist
Acceptance criteria:
- PASS: `A` computed via `CRM_NLEquation_AnalyticalJac` and reshaped from `out_fjac` using the canonical MINPACK column-major conversion helper.
- PASS: `B` produced using raw tip IVP block `x_N._u_zc` only, with no eliminated Jacobians and no `pseudoInverse()`.
- PASS: FD check for `B` passes at one operating point with printed tolerances (see Test Evidence).

Audit checklist:
- PASS: `A` comes from `CRM_NLEquation_AnalyticalJac` (no numerical differencing).
- PASS: `B` uses PR-1 raw Jacobian `x_N._u_zc` and mirrors `K_last` vs identity logic from `CRM_NLEquation_AnalyticalJac`.
- PASS: FD test perturbs only currents `u` and holds `deltau0` fixed.
- PASS: No contact-mode (`FIXED_TIP`) paths invoked.

## 4) Mathematical / Autodiff Details
- Jacobians introduced:
  - `A = ∂F_scaled/∂x_scaled` from `CRM_NLEquation_AnalyticalJac` (free-tip).
  - `B = ∂F_scaled/∂u` from raw tip IVP block `x_N._u_zc`.
- Linear systems solved: none.
- MINPACK/solver layout handling: `MinpackJacobianToRowMajor(fjac, rows, cols, out_rm)` used for all `out_fjac` conversions. Canonical rule `J_rm[r*cols+c] = fjac[r + c*rows]`.

## 5) Tests & Evidence
Commands run:
- `cmake -S . -B build_v1`
- `cmake --build build_v1 -j`
- `c++ -std=c++17 -I./src -I./numerical -I./main build_v1/pr2_fd_driver.cpp main/CRM_KinematicsTestFunctions.cpp build_v1/libCRMCPPLib.a -o build_v1/pr2_fd_driver`
- `./pr2_fd_driver` (run from `build_v1/` so data paths resolve)

Results:
- FD test output:
  - `B FD max_abs_err: 3.61855e-10`
  - `B FD rel_err: 7.96242e-11`
  - `solver localmin: 0`

## 6) Determinism & Safety Verification
- Confirmed: no `pseudoInverse()`, no SVD/damping/regularization introduced in PR-2 paths.
- Confirmed: no differentiation through solver iterations.
- Confirmed: deterministic recomputation uses only forward-cached inputs (no globals/statics added).

## 7) Go / No-Go Decision
GO to CP-2.
