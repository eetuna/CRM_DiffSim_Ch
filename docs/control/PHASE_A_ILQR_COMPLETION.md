# Phase A iLQR Hardening Completion

## Files changed

- `examples/run_ilqr_circle.py`
- `examples/run_ilqr_lemniscate.py`
- `examples/warmup_crm_extension.py`
- `tests/test_ilqr_full_fidelity_smoke.py`
- `tests/test_ilqr_mpc_smoke.py`
- `python/crm_diffsims/control/ilqr.py`
- `python/crm_diffsims/control/__init__.py`
- `docs/control/PHASE_A_ILQR.md`

## Instrumentation added

- Clamp hit rate tracking inside iLQR forward rollouts.
- Per-iteration iLQR console logs with `clamp_rate` and `sign_flip_rate`.
- Per-MPC-step console logs with clamp/sign-flip rates.
- `tip_index_sanity_check` diagnostic for finiteness and small-perturbation sensitivity of `x_t[24:27]`.
- Non-FAST examples now log mean and final tip error per MPC step and save error-vs-time plots.

## Validation harness updates

- Warmup script: `examples/warmup_crm_extension.py` to pre-build the CRM torch extension.
- FAST mode: `--fast` in both example scripts to reduce horizon/steps/iterations and skip plotting.
- Full-fidelity smoke test (real linearization) in `tests/test_ilqr_full_fidelity_smoke.py`.

## Evidence (console excerpts)

- Warmup:
  - `WARMUP OK`
- Full-fidelity smoke test:
  - `iLQR iter 1: cost=1.295914e+03 clamp_rate=0.00% sign_flip_rate=0.00%`
  - `MPC step 1/1: clamp_rate=0.00% sign_flip_rate=0.00%`
  - `1 passed in 113.38s (0:01:53)`
- Circle full run:
  - `tip_index_sanity_check: tip=[-0.6389213036270502, 40.57857203982636, 83.71590997780315] tip_delta_norm=7.679e-02 finite=True`
  - `iLQR iter 1: cost=4.490563e+03 clamp_rate=0.00% sign_flip_rate=0.00%`
  - Timed out after 3600s before completion.
- Lemniscate full run:
  - `tip_index_sanity_check: tip=[-0.6389213036270502, 40.57857203982636, 83.71590997780315] tip_delta_norm=7.679e-02 finite=True`
  - `iLQR iter 1: cost=nan clamp_rate=0.00% sign_flip_rate=0.00%`
  - Timed out after 3600s before completion; repeated `Coil integration Unbounded!!` warnings.

## Commands run

- `python3 examples/warmup_crm_extension.py`
- `python3 -m pytest tests/test_ilqr_full_fidelity_smoke.py -s`
- `python3 examples/run_ilqr_circle.py`
- `python3 examples/run_ilqr_lemniscate.py`

## GO/NO-GO for Phase B (RL)

NO-GO: full-fidelity example runs did not complete within 1 hour, and the full-fidelity smoke test runtime exceeded the <60s target (113s). Lemniscate run also emitted repeated `Coil integration Unbounded!!` warnings and produced NaN cost early. Recommend investigating solver stability/runtime before proceeding to Phase B.
