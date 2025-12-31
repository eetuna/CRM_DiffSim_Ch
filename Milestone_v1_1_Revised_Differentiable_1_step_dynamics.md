# Milestone v1.1 (Revised): Differentiable 1-step free-space dynamics

## 0) Scope (unchanged)
- Implement a differentiable one-step dynamics function: $x_{t+1} = \text{step}(x_t, u_t, \theta)$ for free-space.
- Reuse the repo’s existing stable dynamics stepping path (the same sequence used by the existing dynamics code/MEX wrapper).
- Any embedded nonlinear solve must be differentiated implicitly; do not differentiate through solver iterations.
- v1.1 supports gradients w.r.t. $u_t$ (actuation currents) and $x_t$ (state) only.
- $L_i$ is fixed (no gradient). $\theta$ gradients are not supported.

## 0.1) Implementation Risk Gates
- **MINPACK Jacobian layout (hard gate):** MINPACK-style Jacobians are Fortran column-major. Require a single canonical conversion helper and forbid ad-hoc indexing/reshaping of `fjac`. Reject PR if any callsite converts `fjac` by hand.
- **Determinism (hard gate):** backward must be able to recompute all required Jacobians/VJPs from cached forward data only. Forbid hidden globals/statics/mutable singletons. Reject PR if any backward result depends on non-cached mutable state.
- **Linear solves (hard gate):** no `pseudoInverse()`, no SVD-based solves, no damping/regularization in v1.1. Adjoint solves must be direct dense solves (LU with pivoting; QR acceptable). Reject PR if violated.

## 1) Dynamics state $x_t$ (precise definition; unchanged public layout)
Let $N_{\text{act}} := \text{NUM_ACT_SET}$ (compile-time; v1.1 runtime asserts $N_{	ext{act}} == 1$ but tensors carry the $N_{	ext{act}}$ dimension).

Per batch element, $x_t$ contains:
- $v_L$ (coil linear velocity, coil frame): $[N_{\text{act}}, 3]$
- $w_L$ (coil angular velocity, coil frame): $[N_{\text{act}}, 3]$
- $m_L$ (coil-frame internal moment): $[N_{\text{act}}, 3]$
- $n_L$ (coil-frame internal force): $[N_{\text{act}}, 3]$
- $p_L$ (coil position, spatial frame): $[N_{\text{act}}, 3]$
- $R_L$ (coil orientation, spatial frame): $[N_{\text{act}}, 3, 3]$ stored as 9 ambient values
- $x_{\text{tip}}$ (catheter tip state): $[15]$ packed as $p(3), R(9 \text{ row-major}), u(3)$.

### PyTorch-facing packed state
- $x_t$: shape $[B, 24 \times N_{\text{act}} + 15]$
- Packing order (fixed):
  - $v_L$ ($3 \times N_{\text{act}}$)
  - $w_L$ ($3 \times N_{\text{act}}$)
  - $m_L$ ($3 \times N_{\text{act}}$)
  - $n_L$ ($3 \times N_{\text{act}}$)
  - $p_L$ ($3 \times N_{\text{act}}$)
  - $R_L$ ($9 \times N_{\text{act}}$, row-major per actuator)
  - $x_{\text{tip}}$ (15)

### Batching
- v1.1 runs per-batch element independently (CPU loop).
- Runtime assert: $x_t.\text{shape}[1] == 24 \times N_{\text{act}} + 15$ and $u_t.\text{shape} == [B, N_{\text{act}}, 3]$.
- Runtime assert: $N_{\text{act}} == 1$ for v1.1 build.

### 1.1) Batching semantics disclaimer
- Each batch element is processed independently in forward and backward.
- No coupled, vectorized, or shared-state solves across batch elements are allowed in v1.1.

## 2) Forward API (PyTorch-facing signature; determinism; no randomness)
### Signature
- $x_{tp1} = \text{crm_step}(x_t, u_t, L_i, \text{cfg_dyn}) \rightarrow x_{t+1}$

### Inputs
- $x_t$: $[B, 24 \times N_{\text{act}} + 15]$, dtype float64, requires_grad=True/False as needed.
- $u_t$: $[B, N_{\text{act}}, 3]$, dtype float64, requires_grad=True.
- $L_i$: $[B]$ (or scalar), dtype float64, requires_grad=False.
- $\text{cfg_dyn}$: opaque handle storing all fixed parameters (catheter model/config + dynamics parameters + solver tolerances).

### Output
- $x_{t+1}$: $[B, 24 \times N_{\text{act}} + 15]$ (same packing).

### Determinism
- No randomness.
- Warm-start only via explicit content of $x_t$ (not via globals/statics).
- Must cache and return solver exit status and residual norm (see Section 5).

## 3) Explicit $\eta$ identification (MANDATORY; removes ambiguity)
### 3.1 Requirement
Do not assume $\eta := [m_L, n_L]$ as a mathematical object. Instead:
- Define $\eta$ exactly as the nonlinear decision variable vector passed to the dynamics root-solver routine invoked by `DynamicsBVP` (the `x[]` array provided to `TrustRegionDogleg_dyn`).
- The v1.1 implementation MUST document:
  - $\text{dim}(\eta)$ (expected `NUM_DYN_RESIDUAL` but must be verified)
  - ordering and meaning of each block in $\eta$
  - scaling applied to $\eta$ (e.g., `IVALUE_SCALE_M`, `IVALUE_SCALE_N` style scaling)

### 3.2 Contract after identification
Once the code-level $\eta$ is identified, all Jacobians in v1.1 are defined w.r.t. that $\eta$ exactly (including scaling):
- $A_{\text{step}} := \partial G / \partial \eta$
- $B_{\text{step}} := \partial G / \partial u_t$
- $C_{\text{step}} := \partial G / \partial x_t$
$x_t$ continues to store $m_L, n_L$ as part of the state and as warm-start, but $\eta$ is whatever the solver actually uses.

## 4) Precise definition of residual $G$ (MANDATORY; exact evaluation path)
### 4.1 Requirement
Define $G$ as the exact residual vector evaluated inside the nonlinear solver used in `DynamicsBVP`, with no iteration dependence.
- $G(\eta, x_t, u_t, L_i, \text{cfg_dyn}) \rightarrow r$
- $r$ must match the solver’s residual buffer passed to `TrustRegionDogleg_dyn`.

### 4.2 Evaluation path naming (must match repo call chain)
The plan REQUIRES implementing a callable residual evaluation that follows the same path as forward solve:
- `DynamicsBVP(...)` $\rightarrow$ `TrustRegionDogleg_dyn(...)` $\rightarrow$ calls the solver residual callback (the dynamics NLE function).
- v1.1 requires $G$ be callable directly at a converged $\eta^*$ (no iterations) to compute:
  - residual norm
  - Jacobians $A_{\text{step}}, B_{\text{step}}, C_{\text{step}}$
  - VJP products for the output map (Section 6)

## 5) Determinism & non-convergence policy (MANDATORY)
### 5.1 Forward caching requirements
Forward MUST cache for each batch element:
- $\eta^*$ (the converged decision variable in solver coordinates)
- `solver_exit_flag` (the solver info / status code used by the dynamics solver)
- $\text{residual_norm} = \|G(\eta^*, \dots)\|_2$ (computed consistently)
- any auxiliary values produced as side-outputs of residual evaluation that are required by the forward step mapping (e.g., outputs analogous to `out_u0`, `tau` if the dynamics residual callback produces them)

### 5.2 Backward behavior on non-convergence (explicit)
If `solver_exit_flag` indicates non-convergence OR `residual_norm` exceeds a fixed threshold:
- Default behavior (v1.1): raise an error in backward (and forward may optionally raise as well).
- This policy is required to keep gradcheck deterministic and to avoid silently wrong gradients.
- (If a permissive “zero-grad + warning” mode is desired later, it must be a separate explicit flag in `cfg_dyn`; not part of v1.1 default.)

> **Strict non-convergence policy (hard gate):**
> - If `solver_exit_flag != success` OR `residual_norm > threshold`, backward MUST raise an error (no silent zero-grad).
> - Log (at minimum): `solver_exit_flag`, `residual_norm`, batch index, and time-step index. If iteration count / function-eval count is available from the underlying solver wrapper, log that too.

### 5.3) State overwrite semantics
- $x_{t+1}$ MUST store the converged implicit variable $\eta^*$ (mapped back into public state slots as applicable).
- The initial guess from $x_t$ must never be propagated forward.

## 6) Rotation handling policy (MANDATORY)
$x_t$ contains 9D ambient rotation entries for $R_L$ (and also 9D rotation inside $x_{\text{tip}}$).

Before calling any solver/integrator in forward:
### Assert-orthonormal policy (v1.1 default):
- Check $\|R^T R - I\|_F < \text{tol}_R$ and $\det(R) > 0$ for each stored rotation.
- If violated: treat as non-convergence and follow Section 5.2 behavior.

### Gradient policy:
- Gradients are in ambient $\mathbb{R}^9$ for rotation entries in v1.1 (no Lie tangent mapping, no projection Jacobians).

### 6.1) No hidden projection / stabilization rule
- Aside from the explicit rotation orthonormality check already defined, no projection, clipping, normalization, or stabilization is permitted in forward or backward.
- Any such behavior must be explicit and gated behind future config flags (not part of v1.1).

## 7) Backward specification (implicit differentiation; Jacobian layout contract; VJP products)
### 7.1 Step structure
Forward step is modeled as:
1. Solve implicit nonlinear system for $\eta^*$: $G(\eta^*, x_t, u_t, L_i, \text{cfg_dyn}) = 0$
2. Compute explicit next state: $x_{t+1} = H(\eta^*, x_t, u_t, L_i, \text{cfg_dyn})$

### 7.2 Jacobians required (w.r.t. identified $\eta$)
- $A_{\text{step}} = \partial G / \partial \eta$ (square)
- $B_{\text{step}} = \partial G / \partial u_t$
- $C_{\text{step}} = \partial G / \partial x_t$

### 7.3 Jacobian memory layout contract (MANDATORY)
Any Jacobian returned by MINPACK-style routines (if used anywhere) uses Fortran column-major layout:
- $J(i,j)$ stored at `fjac[i + j*n]`.
v1.1 contracts:
- All internal Jacobians exposed to the PyTorch backward code must be converted to a single canonical layout: row-major dense matrices (or linear-operator applies).
- Transpose conventions in the implicit VJP solve use mathematical transpose of these canonical matrices.
- Canonical internal dense format: row-major with `J_rm[r][c]` stored at `r*cols + c`.
- Canonical row-major reshape rule (required): for a Jacobian with `rows=m`, `cols=n`, the canonical row-major matrix `J_rm` must satisfy `J_rm[r][c] = fjac[r + c*m]`.
- Worked example (2×2):
  - If the mathematical Jacobian is `J = [[a, b], [c, d]]`, then MINPACK `fjac` is `[a, c, b, d]` (column-major).
  - After conversion, canonical row-major `J_rm` must be `[a, b, c, d]`.

### 7.4 Prefer VJP products for $H$ (MANDATORY)
Do not form full $H_x, H_u, H_{\eta}$ matrices by default.
Instead, implement VJP product functions:
- $H_{\eta}^T g = (\partial H / \partial \eta)^T g$
- $H_u^T g = (\partial H / \partial u_t)^T g$
- $H_x^T g = (\partial H / \partial x_t)^T g$
These must be computable at the converged $(\eta^*, x_t, u_t)$ without iterations.

### 7.7 Temporary v1.1 FD bridge (PR-7/CP-7/PR-8 only)
As a temporary bridge for v1.1, FD-based Jacobians/VJP products are allowed for PR-7/CP-7/PR-8 **only**:
- Allowed: `fdjac1_dyn` for $A_{\text{step}}$ and central differences for $B_{\text{step}}, C_{\text{step}}$, and $H_*^T g$.
- Allowed (bridge-only): direct FD of `crm_step_forward` to approximate gradients w.r.t. $x_t$ and $u_t$ as a black-box (no differentiation through solver iterations).
- Not allowed: differentiation through solver iterations, pseudo-inverse/SVD/damping/regularization, or any hidden state.
- Determinism and the strict non-convergence policy remain mandatory and unchanged.
- Any MINPACK layout conversions must use the canonical `MinpackJacobianToRowMajor` helper.
- **Performance/accuracy caveat:** FD-based Jacobians/VJPs are slower and less accurate; results may be sensitive to step size. They are accepted for v1.1 only and must be replaced in v1.2.

### 7.5 Exact implicit VJP solve (transpose conventions)
Given upstream gradient $g = \partial L / \partial x_{t+1}$ (shape $[24 \times N_{\text{act}} + 15]$):
1. Compute RHS using VJP product: $\text{rhs} = H_{\eta}^T g$
2. Solve adjoint: $A_{\text{step}}^T \lambda = \text{rhs}$
3. Compute gradients:
   - $g_u = H_u^T g - B_{\text{step}}^T \lambda$
   - $g_x = H_x^T g - C_{\text{step}}^T \lambda$
No differentiation through solver iterations is allowed.

### 7.6) Adjoint linear solve policy
- The adjoint system $A_{\text{step}}^T \lambda = \text{rhs}$ must be solved using a direct dense solver (LU with pivoting).
- No pseudo-inverse, SVD, damping, or regularization is allowed in v1.1.
- Failure or numerical singularity is treated as non-convergence and follows the existing backward error policy.

## 8) Jacobian sources (file:function) and what must be added (no FD in backward)
For each required Jacobian/VJP block, the plan requires a concrete source:
- $G$ evaluation: MUST be callable through the same residual callback used inside the nonlinear solver invoked by `DynamicsBVP`.
  - Source: existing dynamics residual function used by the solver (repo lookup required; see Section 10).
- $A_{\text{step}}, B_{\text{step}}, C_{\text{step}}$:
  - v1.1 bridge: FD is allowed in backward for PR-7/CP-7/PR-8 (see Section 7.7).
  - Analytic implementation is mandatory for v1.2 (see Section 11).
- $H_{\eta}^T g, H_u^T g, H_x^T g$:
  - MUST be implemented as VJP products (default).
  - v1.1 bridge: FD-based VJP products are allowed for PR-7/CP-7/PR-8 (see Section 7.7).
  - No full-matrix construction required unless trivially available.

## 9) Validation gate for v1.1 (FD tests + rollout smoke test; explicit tolerances)
### 9.1 FD tests for $\partial x_{t+1} / \partial u_t$ ($\\geq 3$ operating points)
- Use float64, B=1, $N_{\text{act}}=1, L_i$ fixed (e.g., 80).
- Define 3 action operating points:
  - $[0.0, 0.0, 0.05]$
  - $[0.02, -0.01, 0.08]$
  - $[-0.03, 0.015, 0.06]$
- Use central differences on $u_t$ with $h=1\text{e-}6$ (fallback $1\text{e-}5$).
- Compare gradients for a scalar loss $L$ that depends on a fixed subset of $x_{t+1}$ (e.g., tip position entries and coil position entries).
- Tolerances:
  - $\text{max_abs_err} \leq 1\text{e-}4$
  - $\text{rel_err} \leq 1\text{e-}3$

### 9.2 Rollout smoke test (10–50 steps)
Run $T=10$ steps with bounded $u_t$ sequence.
- Hard checks:
  - No NaNs in $x_t$.
  - All steps converge (`solver_exit_flag` indicates success).
  - `residual_norm` stays below threshold.
- Log:
  - runtime per step
  - max residual norm

### 9.3 Commands
- Build: `cmake -S . -B build && cmake --build build -j`
- Tests:
  - `test_crm_step_fd.py`
  - `test_crm_step_rollout_smoke.py`
  - `test_crm_step_gradcheck.py`

## 10) Deliverables + checkpoints (PR-style; unchanged structure)
### PR-5: v1.1 Forward step() wrapper (existing dynamics stepping path)
- **Goal:** Implement deterministic `step_forward` using the repo’s existing dynamics pipeline.
- **Tasks:** Implement pack/unpack; call the existing dynamics step path; cache $\eta^*$, exit flag, residual norm.
- **Acceptance:** Deterministic repeatability; $N_{\text{act}}==1$ enforced; forward returns $[B, 24 \times N_{\text{act}}+15]$.
- **Checkpoint:** CP-5 (Forward step stable + cached $\eta^*$/flags)

### PR-6: v1.1 $\eta$ identification + residual $G$ callable + Jacobian layout contract
- **Goal:** Lock the exact $\eta$ definition and implement $G(\eta, \dots)$ callable at $\eta^*$ with a documented Jacobian layout contract.
- **Tasks:** Inspect `DynamicsBVP`/solver call chain; document $\eta$ dimension/order/scaling; implement $G$ evaluation at $\eta^*$; define canonical row-major Jacobian representation.
- **Acceptance:** $\eta$ contract written and validated by printing $\eta^*$ and residual norm; layout conversions documented.
- **Checkpoint:** CP-6 ( $\eta$ + $G$ contract locked)

#### $\eta$ (eta) contract (code-locked; v1.1)
Define $\eta$ as **exactly** the decision vector passed to `TrustRegionDogleg_dyn(n, x, ...)` inside `DynamicsBVP(...)` (i.e., the solver’s `double x[]` in solver coordinates).
- **Dimension:** $\dim(\eta) = \text{NUM\_DYN\_RESIDUAL} = 6 \times \text{NUM\_ACT\_SET}$ (for v1.1, runtime asserts `NUM_ACT_SET == 1`, so $\dim(\eta)=6$).
- **Ordering (per actuator set `j`):** the 6-vector block is `[m_L_scaled (3), n_L_scaled (3)]`:
  - `eta[6*j + 0..2]` = `m_L_scaled[j][0..2]`
  - `eta[6*j + 3..5]` = `n_L_scaled[j][0..2]`
- **Scaling rules (must match forward solver code path):**
  - Packing in `DynamicsBVP(...)` uses:
    - `m_L_scaled = (1 / IVALUE_SCALE_M) * m_L_initialguess`
    - `n_L_scaled = (1 / IVALUE_SCALE_N) * n_L_initialguess`
  - Unpacking in `DYNNLEquation(...)` uses:
    - `m_L = IVALUE_SCALE_M * m_L_scaled`
    - `n_L = IVALUE_SCALE_N * n_L_scaled`
- **Pack/unpack definitions (by function name / file):**
  - Pack to solver coordinates: `DynamicsBVP(...)` in `src/CoilDynamics_Defs.cpp` (the `initialguessscaled[...]` / `x[...]` loops).
  - Unpack from solver coordinates: `DYNNLEquation(...)` in `src/CoilDynamics_Defs.cpp` (the `m_L[j][i] = IVALUE_SCALE_M * in_x[...]` and `n_L[j][i] = IVALUE_SCALE_N * in_x[...]` loops).

### PR-7: v1.1 Analytic Jacobians for $A_{\text{step}}, B_{\text{step}}, C_{\text{step}}$ and VJP products for $H$
- **Goal:** Provide analytic $A_{\text{step}}, B_{\text{step}}, C_{\text{step}}$ and VJP products $H_*^T g$ (no FD in backward).
- **Tasks:** Implement analytic Jacobians consistent with $\eta$ scaling and residual definition; implement $H_{\eta}^T g, H_u^T g, H_x^T g$.
- **Acceptance:** Jacobian/VJP functions evaluate at $\eta^*$ with consistent dimensions and no iteration dependence.
- **Checkpoint:** CP-7 (Jacobians/VJPs available)

### PR-8: PyTorch autograd for step() + validation suite
- **Goal:** Integrate implicit VJP backward in PyTorch and pass FD/gradcheck + rollout smoke tests.
- **Tasks:** Bind `step_forward`; implement backward solve $A_{\text{step}}^T \lambda = H_{\eta}^T g$; compute $g_u, g_x$; enforce non-convergence policy; add tests.
- **Acceptance:** gradcheck passes in stable regime; FD tests pass tolerances; rollout smoke test passes.
- **Checkpoint:** CP-8 (v1.1 complete)

## 11) Milestone v1.2: Replace FD Jacobians/VJPs with analytic implementations
### Goal
Replace the v1.1 FD bridge with analytic $A_{\text{step}}, B_{\text{step}}, C_{\text{step}}$ and analytic $H_*^T g$ products.

### Acceptance criteria
- Analytic $A_{\text{step}}$ from the dynamics residual Jacobian callback used by `TrustRegionDogleg_dyn` (solver coordinates).
- Analytic $B_{\text{step}}$ and $C_{\text{step}}$ assembled via IVP analytic Jacobian blocks and chain rule.
- Analytic $H_{\eta}^T g$, $H_u^T g$, $H_x^T g$ via analytic IVP Jacobians and chain rule (no FD).
- No differentiation through solver iterations; strict non-convergence policy unchanged.
- Canonical MINPACK column-major to row-major conversion helper used wherever applicable.

### Tests
- FD comparison with tighter tolerances than v1.1 (target: `max_abs_err <= 1e-5`, `rel_err <= 1e-4`).
- Gradcheck with tightened tolerances (target: `rtol <= 1e-4`, `atol <= 1e-6`) in stable regimes.
- Runtime benchmark: report per-step runtime with analytic Jacobians vs FD bridge.

## Minimal Python training-loop demo (pseudo-code; batched; gradients to $u_t$)
```python
import torch

B, T, N_act = 8, 20, 1
cfg_dyn = load_cfg_dyn(...)
Li = torch.full((B,), 80.0, dtype=torch.float64, requires_grad=False)

x = x0.clone().detach()                  # [B, 24*N_act + 15]
u = torch.zeros((T, B, N_act, 3), dtype=torch.float64, requires_grad=True)

loss = 0.0
for t in range(T):
    x = crm_step(x, u[t], Li, cfg_dyn)   # differentiable
    tip_p = x[:, 24*N_act : 24*N_act+3]  # tip position inside x_tip
    loss = loss + (tip_p**2).sum(dim=1).mean()

loss.backward()
# u.grad is available for MPC/iLQR/model-based RL updates
```

## Remaining repo lookups required
1. **Confirm solver status / info semantics:**
   - Inspect `minpack_DYN_Defs.cpp` implementation of `TrustRegionDogleg_dyn` to confirm how `info` is set and how “success” vs “non-convergence” is determined for v1.1 `solver_exit_flag`.
2. **Confirm the forward stepping path equivalence to existing stable pipeline:**
   - Inspect `CRMDYN_c_mex.cpp` and/or `CRMDYN_c.cpp` to confirm the exact order of calls and state assembly used in the current step function.
3. **Confirm rotation validity expectations:**
   - Inspect where `R_pre`/`R_L` are produced and whether they are expected orthonormal by construction in the existing integrator path.

## Doc change summary
- Added explicit implementation risk gates (layout conversion, determinism, linear-solve restrictions) with reject-if-violated criteria.
- Made the MINPACK Jacobian layout contract explicit with a canonical reshape rule and worked example.
- Added a code-locked $\eta$ contract (dimension/order/scaling) tied to the actual `DynamicsBVP`/`TrustRegionDogleg_dyn` decision vector and identified the pack/unpack locations by function name/file.
- Added a strict non-convergence policy callout with required logging for debugging.
