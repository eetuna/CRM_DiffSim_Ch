# PR-9 Completion

## Files Changed
- python/crm_diffsims/__init__.py
- python/crm_diffsims/dynamics/__init__.py
- python/crm_diffsims/dynamics/step.py
- pyproject.toml
- examples/run_step_rollout_and_backward.py
- tests/test_crm_step_fd.py
- tests/test_crm_step_gradcheck.py
- tests/test_crm_step_rollout_smoke.py
- docs/autodiff_audits/PR-9_completion.md
- docs/autodiff_audits/CP-9_completion.md

## Files Removed
- tests/crm_step_common.py
- tests/crm_step_autograd.py

## What Was Moved
- Dynamic extension loader, config builder, state builders, and FD autograd bridge moved from:
  - tests/crm_step_common.py
  - tests/crm_step_autograd.py
- New stable API at python/crm_diffsims/dynamics/step.py with public entry point `crm_step`.

## Behavior Confirmation
- FD forward/backward behavior preserved by moving code verbatim into the package and keeping solver exit/residual checks unchanged.
- No changes to convergence criteria, residual thresholds, or FD epsilon.

## Compliance Checklist
- [x] Public API: `from crm_diffsims.dynamics import crm_step`
- [x] No differentiation through solver iterations (FD only)
- [x] No pseudoInverse/SVD/damping/regularization added
- [x] Deterministic forward/backward
- [x] Strict nonconvergence policy preserved (error on failure)
