# PR-18 Completion

## Files Modified
- tests/test_stress_replay_npz.py
- examples/run_workspace_stress_report.py
- python/crm_diffsims/dynamics/step.py
- docs/autodiff_audits/PR-18_completion.md
- docs/autodiff_audits/CP-18_completion.md

## Equilibrium Seed x0 Construction
- Calls v1.0 equilibrium forward with u0 = currents[0] and Li from dataset.
- Builds x0 as a packed x_t tensor with:
  - v_L_pre, w_L_pre = 0
  - mL_guess, nL_guess = 0
  - p_pre = equilibrium p_tip
  - R_pre = identity (fallback; equilibrium forward does not expose orientation)
  - xf_pre = [p_tip, I_3, 0] (position + identity rotation + zero curvature)
- Deterministic: no randomness, no hidden state.

## Rotation Fallback
- Equilibrium forward exposes only p_tip, so R_pre and xf_pre rotation use identity as a documented fallback.

## Classification Changes
- Added classification:
  - expected_initial_state_mismatch when large tip_error at indices 0–5 with low residual and equilibrium_seed reduces error sharply.
  - expected_nonconvergence when large error coincides with solver failure/high residual.
- Removed “inconsistency” outliers for early indices when equilibrium_seed explains the mismatch.

## Notes
- Added separate equilibrium extension loader (crm_equilibrium_ext) to avoid linking multiple PYBIND modules into a single extension.
- Stress report CSV now includes init_mode for mode-by-mode comparison.
