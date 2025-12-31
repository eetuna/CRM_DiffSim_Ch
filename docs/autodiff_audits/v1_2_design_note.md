# v1.2 Design Note: Analytic Jacobians/VJPs for Dynamics Residual

## Scope
This note documents the exact dynamics residual G(eta, x_t, u_t, Li) as implemented in the v1.1 DYNNLEquation path, identifies dependencies, and specifies analytic Jacobian/VJP blocks to implement for v1.2. No behavior changes are proposed here.

## Residual Definition (DYNNLEquation)
The residual vector is G(eta, x_t, u_t, Li) = out_y in `DYNNLEquation` (`src/CoilDynamics_Defs.cpp:427`).

Inputs:
- eta = stacked [m_L_scaled, n_L_scaled] for each actuator set.
- x_t = packed state (v_L_pre, w_L_pre, mL_guess, nL_guess, p_pre, R_pre, xf_pre) in `crm_step_forward` and `crm_dyn_residual` (`src/CRM_TorchDynamics.cpp:226` and `src/CRM_TorchDynamics.cpp:403`).
- u_t = actuator currents used to construct CRMShootingMethodParams (`CRMDYNConstructShootingMethodParamSet`, `src/CRM_TorchDynamics.cpp:300`).
- Li = inserted length used in CRMShootingMethodParams (`src/CRM_TorchDynamics.cpp:300`).

Scaling:
- In DYNNLEquation, eta is scaled into physical m_L, n_L using:
  - m_L = IVALUE_SCALE_M * eta[0:3]
  - n_L = IVALUE_SCALE_N * eta[3:6]
  - See `src/CoilDynamics_Defs.cpp:434` and `src/CRMDYN.hpp:37`.
- Residual entries are scaled with:
  - RESIDUAL_SCALE_P for position residuals
  - RESIDUAL_SCALE_R for rotation residuals
  - See `src/CoilDynamics_Defs.cpp:692` and `src/CRMDYN.hpp:41`.

Residual construction (per actuator set, NUM_ACT_SET == 1 for v1.1):
- The residual is a 6-vector per actuator: [position_residual(3), rotation_residual(3)].
- The algorithm integrates from distal to proximal across segments (flexible and rigid) and computes:
  - p_f, R_f from CRMFlexible_IVP_Back (backward integration over flexible segment)
  - p_L, R_L from coil dynamics integration (CoilDynamics) and rigid segment offset.
- For intermediate coil-to-flex transitions (actno_mn = actno + 1):
  - position residual: p_f - p_L
  - rotation residual: sqrt(vNormSq of column-wise differences) where
    v1 = R_f(:,0) - R_L(:,0)
    v2 = R_f(:,1) - R_L(:,1)
    v3 = R_f(:,2) - R_L(:,2)
    residual_rot = [sqrt(||v1||^2), sqrt(||v2||^2), sqrt(||v3||^2)]
  - See `src/CoilDynamics_Defs.cpp:630-679`.
- For the final (proximal) residual (actno = 0):
  - position residual: p_f - p_d, where p_d is root pose position (Params.xi)
  - rotation residual: same column-wise norm residual between R_f and R_d (Params.xi)
  - See `src/CoilDynamics_Defs.cpp:688-680` and `src/CoilDynamics_Defs.cpp:688-683`.
- Output assembly:
  - out_y[0:2] = RESIDUAL_SCALE_P * position residual
  - out_y[3:5] = RESIDUAL_SCALE_R * rotation residual
  - See `src/CoilDynamics_Defs.cpp:692-697`.

Exact dependencies inside DYNNLEquation:
- eta (m_L_scaled, n_L_scaled):
  - Directly used to form m_L and n_L (scaled) (`src/CoilDynamics_Defs.cpp:434-439`).
  - m_L drives flexible segment curvature u_L (via Kinv, ustar) and net moment net_mL.
  - n_L drives net_nL used by CoilDynamics.
- x_t:
  - v_L_pre, w_L_pre, p_pre, R_pre are loaded into x_coil initial state (`src/CoilDynamics_Defs.cpp:455-472`).
  - xf_pre is used as tip pose for the last flexible segment (Params.xf) (`src/CRM_TorchDynamics.cpp:149`).
  - mL_initialguess, nL_initialguess are passed into CRMDYNSolverIVP_Prep, but do not enter DYNNLEquation directly; they affect params and backward IVP calls.
- u_t:
  - enters CRMShootingMethodParams via CRMDYNConstructShootingMethodParamSet; affects magnetic moment, coil alignment, and IVP parameters used by DYNNLEquation and CRMFlexible_IVP_Back/CoilDynamics (`src/CRM_TorchDynamics.cpp:300`, `src/CRMDYN.hpp:62`).
- Li:
  - enters CRMShootingMethodParams (InsertedLength and segment bounds), affecting segment geometry and IVP setup used by DYNNLEquation (`src/CRM_TorchDynamics.cpp:300`).

## Jacobian Blocks to Implement
We target analytic Jacobians of G with respect to eta, u_t, and x_t.

Definitions (per batch b):
- A_step = dG/deta, shape [NUM_DYN_RESIDUAL, NUM_DYN_RESIDUAL].
- B_step = dG/du_t, shape [NUM_DYN_RESIDUAL, 3 * NUM_ACT_SET].
- C_step = dG/dx_t, shape [NUM_DYN_RESIDUAL, 24 * NUM_ACT_SET + 15].
- G is the residual from `DYNNLEquation` (`src/CoilDynamics_Defs.cpp:427`).
- For v1.1, NUM_ACT_SET = 1, so NUM_DYN_RESIDUAL = 6.

Transpose convention:
- Jacobian blocks are row-major with residual rows and input columns:
  - A_step[i,j] = dG_i / d(eta_j)
  - B_step[i,j] = dG_i / d(u_t_j)
  - C_step[i,j] = dG_i / d(x_t_j)
- This matches existing FD outputs in `crm_dyn_jacobians` (`src/CRM_TorchDynamics.cpp:520`).

IVALUE_SCALE_M/N involvement:
- eta is scaled before entering physics:
  - d(m_L)/d(eta_m) = IVALUE_SCALE_M * I_3
  - d(n_L)/d(eta_n) = IVALUE_SCALE_N * I_3
- These scaling factors must be included in analytic A_step.
- Scaling occurs in `DYNNLEquation` and in `EvaluateDynHOutput` for H (see below).

## VJP Products to Implement
We target analytic products for H, the step mapping used for x_tp1 generation in `EvaluateDynHOutput` (`src/CRM_TorchDynamics.cpp:122`).

Definition:
- H(eta, x_t, u_t, Li) := x_tp1, where H is the mapping implemented by `EvaluateDynHOutput`.
- Required VJPs for a given upstream gradient g = dL/dx_tp1:
  - H_eta^T g, shape [NUM_DYN_RESIDUAL]
  - H_u^T g, shape [3 * NUM_ACT_SET]
  - H_x^T g, shape [24 * NUM_ACT_SET + 15]

Transpose convention:
- H_eta^T g means (dH/deta)^T * g, where g is a column vector with size [24 * NUM_ACT_SET + 15].
- H_u^T g and H_x^T g follow the same convention.
- Matches current FD VJP layout in `crm_dyn_vjp` (`src/CRM_TorchDynamics.cpp:790`).

IVALUE_SCALE_M/N involvement for H:
- In `EvaluateDynHOutput`, eta_scaled is mapped back to physical m_L, n_L via IVALUE_SCALE_M/N (`src/CRM_TorchDynamics.cpp:142`).
- The analytic H_eta^T g must include these scales.

## Residual Term Dependencies (Summary)
Residual entries G are functions of:
- eta:
  - m_L_scaled and n_L_scaled affect u_L, net_mL, net_nL, and thus p_f, R_f and coil states.
- x_t:
  - v_L_pre, w_L_pre, p_pre, R_pre set initial coil states.
  - xf_pre sets distal pose constraints for backward flexible integration.
  - mL_initialguess, nL_initialguess affect CRMDYNSolverIVP_Prep and flexible integration setup.
- u_t:
  - affects magnetic moments and actuation parameters used in CoilDynamics and IVP setup through CRMShootingMethodParams.
- Li:
  - affects segment bounds, insertion length, and IVP setup.

## Source Mapping
Key locations to reference for analytic derivation:
- Residual definition and assembly:
  - `src/CoilDynamics_Defs.cpp:427` (DYNNLEquation)
  - `src/CoilDynamics_Defs.cpp:630` (intermediate position residual)
  - `src/CoilDynamics_Defs.cpp:655` (rotation residuals)
  - `src/CoilDynamics_Defs.cpp:688` (proximal residual)
  - `src/CoilDynamics_Defs.cpp:692` (RESIDUAL_SCALE_P/R)
- Parameter construction and inputs:
  - `src/CRM_TorchDynamics.cpp:300` (CRMDYNConstructShootingMethodParamSet)
  - `src/CRM_TorchDynamics.cpp:403` (crm_dyn_residual input unpacking)
- Step output mapping (H):
  - `src/CRM_TorchDynamics.cpp:122` (EvaluateDynHOutput)
  - `src/CRM_TorchDynamics.cpp:167` (m_L/n_L from eta_scaled)
- Scaling definitions:
  - `src/CRMDYN.hpp:37` (IVALUE_SCALE_M/N)
  - `src/CRMDYN.hpp:41` (RESIDUAL_SCALE_P/R)

## Reuse vs New Derivations
Existing analytic Jacobians (reusable pieces):
- IVP Jacobians for tip kinematics are available in:
  - `CRMSolverIVPJacobian` and `CRMSolverIVPJacobian_TipRaw` (`src/CRM_IVPJacobian.cpp:1`).
- These provide derivatives of tip pose and strain with respect to deltau0 and actuation vectors for the *forward* IVP.

New derivations required for v1.2 dynamics residual:
- DYNNLEquation residual uses backward flexible integration (CRMFlexible_IVP_Back) and CoilDynamics, which currently lack analytic Jacobians.
- A_step (dG/deta), B_step (dG/du_t), C_step (dG/dx_t) must be derived for:
  - CRMFlexible_IVP_Back contributions to p_f, R_f
  - CoilDynamics contributions to out_x_coil and p_L/R_L
  - Composition of residual terms (p_f - p_L, column-difference norms)
- VJP for H (EvaluateDynHOutput) requires derivatives through:
  - EvaluateDynResidual (to obtain u0 and tau) and DYNSolverIVP
  - Mapping from eta_scaled to m_L/n_L with IVALUE_SCALE_M/N
- Existing IVP Jacobians may be used to validate or reuse sub-blocks for forward IVP terms, but there is no direct analytic Jacobian available for the backward IVP and coil dynamics path used by DYNNLEquation.

## Dimensions and Conventions
For NUM_ACT_SET = 1 (v1.1):
- eta: [6] = [m_L_scaled(3), n_L_scaled(3)]
- x_t: [24 * 1 + 15] = [39]
- u_t: [3]
- G: [6]
- H output (x_tp1): [39]

General shapes:
- A_step: [6 * NUM_ACT_SET, 6 * NUM_ACT_SET]
- B_step: [6 * NUM_ACT_SET, 3 * NUM_ACT_SET]
- C_step: [6 * NUM_ACT_SET, 24 * NUM_ACT_SET + 15]
- H_eta^T g: [6 * NUM_ACT_SET]
- H_u^T g: [3 * NUM_ACT_SET]
- H_x^T g: [24 * NUM_ACT_SET + 15]

## v1.2 Implementation Notes (Planning Only)
- Implement analytic residual Jacobians at the DYNNLEquation layer (not inside the nonlinear solver) to preserve the no-differentiation-through-iterations constraint.
- Keep the same residual scaling, convergence policy, and deterministic behavior as v1.1.
