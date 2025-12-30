Title:
“CRM v1.0 Differentiable Equilibrium – Execution Plan (Frozen)”

Implementation Risk Gates (apply to ALL PRs)
- MINPACK Jacobian layout: MINPACK-style `out_fjac` is Fortran column-major; require a single canonical conversion helper and forbid ad-hoc indexing. Reject PR if any callsite reshapes `out_fjac` by hand.
- Determinism: backward must be able to recompute all Jacobians from cached forward data only; forbid hidden globals/statics/mutable singletons in forward/backward paths. Reject PR if any backward result depends on non-cached mutable state.
- Linear solves: no `pseudoInverse()`, no SVD-based solves, no damping/regularization in v1.0. Solve all adjoint systems with a direct dense solver (LU with pivoting; QR acceptable). Reject PR if violated.

Jacobian Layout Contract (MINPACK `out_fjac` → canonical row-major)
- Canonical internal format for all Jacobians used in backward: row-major dense matrix `J_rm` with `J_rm[r][c]` stored at `r*cols + c`.
- MINPACK/Fortran storage: `fjac` stores `J` column-major, i.e., `fjac[r + c*rows] = J(r,c)`.
- Required conversion rule: `J_rm[r][c] = fjac[r + c*rows]`.
- Worked example (2×2):
  - If the mathematical Jacobian is `J = [[a, b], [c, d]]`, then MINPACK `fjac` is `[a, c, b, d]` (column-major).
  - After conversion, canonical row-major `J_rm` must be `[a, b, c, d]`.

PR-1: Raw IVP Jacobian Getter (TIP)
1) Goal (1 sentence)
- Expose raw, pre-elimination IVP Jacobian blocks at the tip for `u0` and `zc` (currents) without using `pseudoInverse()`.

2) Exact tasks (bullet list, imperative)
- Add an internal C++ function that runs `CRMSolverIVP_CoreWithJacobian(...)` and returns the raw tip Jacobian blocks: `_p_u0`, `_p_zc`, `_u_u0`, `_u_zc`, `_ws_u0`, `_ws_zc`.
- Ensure the function does not compute or return any eliminated BVP Jacobians (`JBVP_*`) and never calls `pseudoInverse()`.
- Confirm the returned matrices use row-major mapping consistent with `CRM_StateVector_Definitions.hpp`.
- Add a minimal C++ test harness that prints/returns these matrices for a fixed parameter set and input currents (no numeric assertions yet).

3) Files to touch (exact paths)
- `src/CRM_IVPJacobian.hpp`
- `src/CRM_IVPJacobian.cpp`
- `main/CRM_KinematicsTestFunctions.cpp` (or add a new minimal test file under `main/` if preferred)

4) Acceptance criteria (checkable)
- The new getter returns the six raw 3×3 blocks at the tip, and each block matches the confirmed member paths: `x_N._p_u0`, `x_N._p_zc`, `x_N._u_u0`, `x_N._u_zc`, `x_N._ws_u0`, `x_N._ws_zc`.
- No callsites of `pseudoInverse()` are introduced or used in the new getter path.
- For `NUM_ACT_SET=1`, the getter returns `zc` as a 3×3 matrix (no `zl` columns included).

5) Audit checklist (YES/NO items)
- YES/NO: Getter uses `CRMSolverIVP_CoreWithJacobian` and reads `x_N._*` members directly.
- YES/NO: Getter never calls `CRMSolverIVPJacobian(...)` (the function that builds eliminated Jacobians).
- YES/NO: Getter path contains zero uses of `pseudoInverse()` or `completeOrthogonalDecomposition()`.
- YES/NO: Returned Jacobians correspond to `zc` for currents (not `zl`).

PR-2: Residual Jacobians A=F_x_scaled, B=F_u_scaled + FD test
1) Goal (1 sentence)
- Compute the scaled equilibrium residual Jacobians `A = ∂F_scaled/∂x_scaled` and `B = ∂F_scaled/∂u` at the converged equilibrium, and validate `B` by finite differences.

2) Exact tasks (bullet list, imperative)
- Trace and lock the exact residual `F_scaled(...)` used in the actual equilibrium solve path:
  - Solver path: `CRMShootingMethodBVP(...)` → `TrustRegionDogleg_GivenJacobian(...)` → `CRM_NLEquation(...)` (free-tip) in `src/CRM_BVPSolver.cpp`.
  - Residual definition (free-tip): `F_scaled = RESIDUAL_SCALE_M * MomentResidual`, where `MomentResidual` is the `MomentResidual[3]` output of `CRMSolverIVP_Core(...)` (invoked by `CRM_NLEquation(...)`).
  - If multiple residual forms exist in the repo, use ONLY the one actually invoked by `CRMShootingMethodBVP` in the `FREE_TIP` path above.
- Add a function to compute `A` by calling `CRM_NLEquation_AnalyticalJac` (free-tip path only) and reshaping `out_fjac` into a 3×3 matrix with the correct column-major interpretation (minpack convention).
- Add a function to compute `B = ∂F_scaled/∂u` using the raw tip IVP Jacobians from PR-1 and the same residual definition used in the equilibrium solve (free-tip moment residual).
  - Required raw tip block(s) for v1.0 `FREE_TIP` moment residual: `x_N._u_zc` only.
  - Assembly rule must mirror the stiffness-mapping logic used by `CRM_NLEquation_AnalyticalJac(...)` for `A`:
    - If `Params.SegmentTypes[Params.no_segments - 1] == FLEXIBLE`: `B = RESIDUAL_SCALE_M * K_last * J_u_zc_tip`
    - Else (last segment rigid): `B = RESIDUAL_SCALE_M * I * J_u_zc_tip` (i.e., `B = RESIDUAL_SCALE_M * J_u_zc_tip`)
    - `J_u_zc_tip` is the 3×(3*NUM_ACT_SET) row-major block at the tip: `x_N._u_zc` from `AugmentedStateVector<IVPJacobiansFull> x_N`.
    - For v1.0 (`NUM_ACT_SET=1`): `J_u_zc_tip` and `B` are 3×3 and must correspond to `zc` (currents) columns only (no `zl`).
- Explicitly do NOT use `x_N._p_zc` or `x_N._ws_zc` when assembling `B` for v1.0 `FREE_TIP` moment residual.
- Add a finite-difference test that perturbs `u` (currents) at fixed `x_scaled` (hold `deltau0` fixed) and compares FD to the computed `B`.
- Ensure the FD test runs with a fixed parameter set and fixed `Li`, and reports max-abs and relative errors.

3) Files to touch (exact paths)
- `src/CRM_BVPSolver.cpp`
- `src/CRM_BVPIVP_APIDeclarations.hpp` (only if a new internal API signature is declared here)
- `main/CRM_KinematicsTestFunctions.cpp` (or a new dedicated test file under `main/`)

4) Acceptance criteria (checkable)
- `A` is computed via `CRM_NLEquation_AnalyticalJac` with correct reshaping of `out_fjac` (Fortran-style column-major) into a C++ 3×3.
- `B` is produced without any use of eliminated Jacobians and without `pseudoInverse()`.
- FD check for `B` passes at at least one operating point with explicit printed tolerances met (as agreed in the frozen checklist).

5) Audit checklist (YES/NO items)
- YES/NO: `A` comes from `CRM_NLEquation_AnalyticalJac` and not from numerical differencing.
- YES/NO: `B` uses PR-1 raw Jacobian `x_N._u_zc` (and the same last-segment `K_last` vs identity logic as `CRM_NLEquation_AnalyticalJac`) and not eliminated Jacobians.
- YES/NO: FD test perturbs only currents `u` and holds `deltau0` fixed.
- YES/NO: No contact-mode (`FIXED_TIP`) paths are invoked.

PR-3: PyTorch Autograd Backward (TIP_POSITION_ONLY)
1) Goal (1 sentence)
- Implement a PyTorch custom autograd function for free-space equilibrium that returns tip position only and backpropagates gradients w.r.t. currents via the frozen implicit VJP solve.

2) Exact tasks (bullet list, imperative)
- Add a C++ extension entry point that runs equilibrium forward (solve `deltau0` via `CRMShootingMethodBVP`, then IVP for `p_tip`) and returns `p_tip` for each batch element.
- Save (cache) the minimal tensors/aux data required by the frozen checklist: `u`, `Li`, converged `deltau0`, `ftip`, and enough parameter/config state to recompute `A`, `B`, `y_x`, `y_u_direct` deterministically.
- Forward caches for backward (required; reject PR if incomplete):
  - | Cache name | dtype / shape | Why needed |
    |---|---:|---|
    | `u` | float64 `[B, N_act, 3]` | Recompute `B` and `y_u_direct`; return `∂L/∂u` |
    | `Li` | float64 `[B]` (or scalar) | Deterministic forward recomputation for `A/B/y_*` (no grad) |
    | `deltau0` | float64 `[B, 3]` | Defines the converged equilibrium state `x_scaled`; needed for `A/B/y_*` |
    | `ftip` | float64 `[B, 3]` | Required to rerun IVP/Jacobian calls consistently (free-tip uses fixed `TipForce`, but cache the effective `ftip` actually used) |
    | `localmin` / solver status | int32 `[B]` | Determinism + backward policy (error on non-convergence) |
    | `cfg_equil` (immutable) | opaque | Must include any solver knobs/scales that affect Jacobians (e.g., `RESIDUAL_SCALE_M`, `IVALUE_SCALE_DU`, solver tolerance), and all physical parameters/config used in IVP/BVP |
- Implement backward for `TIP_POSITION_ONLY` using the frozen implicit VJP:
  - Solve `A^T λ = y_x^T g_p`
  - Return `g_u = y_u_direct^T g_p - B^T λ`
- Add a Python test that runs `torch.autograd.gradcheck` for `TIP_POSITION_ONLY` on 3 operating points (fixed `Li`, currents only require grad).

3) Files to touch (exact paths)
- `CMakeLists.txt` (only as needed to build the PyTorch extension)
- `src/` (new binding source file(s) under `src/`, exact name to be chosen during implementation)
- `src/CRM_ForwardKinematics.cpp` (only if a new non-API helper is required; avoid if possible)
- `src/CRM_BVPSolver.cpp` (only if needed for a clean callable forward wrapper)
- `tests/` (new Python test file; exact path to be created during implementation)

4) Acceptance criteria (checkable)
- Forward returns `p_tip` with shape `[B, 3]` from PyTorch.
- Backward returns gradients only for currents `u` and matches `gradcheck` at agreed tolerances for the agreed operating points.
- No eliminated Jacobians and no `pseudoInverse()` are used in backward.
- `Li` is accepted as input but `requires_grad=False` and produces no gradient.

5) Audit checklist (YES/NO items)
- YES/NO: Autograd backward uses only `A`, `B`, `y_x`, `y_u_direct` and the implicit VJP solve (no differentiation through dogleg iterations).
- YES/NO: `TIP_POSITION_ONLY` path does not compute or save any rotation-related tensors.
- YES/NO: Tensor shapes are expressed as `[B, N_act, 3]` for inputs with runtime check `N_act==1`.
- YES/NO: `pseudoInverse()` appears nowhere in the new backward path.

PR-4 (optional): TIP_POSE rotation tangent mapping + gradcheck
1) Goal (1 sentence)
- Extend the autograd function to optionally include tip pose outputs and propagate rotation gradients via the frozen tangent mapping using `ws` Jacobians.

2) Exact tasks (bullet list, imperative)
- Add `readout_mode` support and implement `TIP_POSE` output: `p_tip` and `R_tip`.
- Save `R_tip` in forward for the tangent mapping in backward.
- In backward, convert upstream `g_R` to `g_φ` using the frozen left-trivialized mapping and reuse the same `A^T` factorization/solve with an additional RHS for the rotation tangent component.
- Compute `y_x` and `y_u_direct` for the rotation tangent using raw `ws` blocks from PR-1 (`_ws_u0` and `_ws_zc`).
- Add `gradcheck` tests for `TIP_POSE` at the agreed operating points and tolerances.

3) Files to touch (exact paths)
- Same PyTorch binding files from PR-3 (under `src/`)
- `tests/` (Python tests)

4) Acceptance criteria (checkable)
- `TIP_POSE` forward returns `(p_tip, R_tip)` with shapes `[B,3]` and `[B,3,3]`.
- `gradcheck` passes for `TIP_POSE` with gradients w.r.t. currents only.
- Rotation gradients use only the `ws` raw Jacobian blocks and the frozen tangent mapping; no log/exp Jacobians are introduced.

5) Audit checklist (YES/NO items)
- YES/NO: Rotation tangent mapping uses left-trivialized convention exactly as frozen.
- YES/NO: Backward adds only extra RHS solves with the same `A^T` (no new solver logic).
- YES/NO: No actuator/coil outputs are added in PR-4.

Doc change summary
- Added explicit implementation risk gates (layout, determinism, linear-solve restrictions) with reject-if-violated criteria.
- Made MINPACK Jacobian layout contract explicit with a worked example to eliminate reshape ambiguity.
- Locked the equilibrium residual definition to the actual `CRMShootingMethodBVP` → `CRM_NLEquation` free-tip code path and made `B` assembly depend explicitly on `x_N._u_zc` only.
- Added a required forward-cache table for autograd backward determinism (what is cached, shapes, and purpose).

STOP. This plan is frozen. Do not implement.
