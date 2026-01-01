# v1.3 PR-1 Completion

## Scope

- Design spec and scaffolding for implicit backward (no math implementation).
- C++ and Python stubs establishing API and cache layout.
- Interface smoke test for shapes/determinism.

## Files changed

- `docs/v1_3/DESIGN.md`
- `src/CRM_TorchDynamics.cpp`
- `python/crm_diffsims/dynamics/step_v1_3.py`
- `tests/test_v1_3_interface_smoke.py`

## Tests run

- `python3 examples/warmup_crm_extension.py`
- `python3 -m pytest tests/test_v1_3_interface_smoke.py -s`

## GO/NO-GO for PR-2

GO: scaffolding and interface tests are in place. Ready to implement implicit backward math and linear solves in PR-2.
