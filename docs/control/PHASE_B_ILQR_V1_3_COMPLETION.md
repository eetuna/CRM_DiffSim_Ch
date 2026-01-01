# Phase B v1.3 iLQR/MPC Completion

## Files changed
- python/crm_diffsims/control/ilqr.py
- examples/run_ilqr_circle.py
- examples/run_ilqr_lemniscate.py
- python/crm_diffsims/dynamics/linearize_v1_3.py
- tests/test_ilqr_mpc_v1_3_circle_smoke.py
- tests/test_ilqr_mpc_v1_3_lemniscate_smoke.py
- scripts/compare_cem_vs_ilqr_v1_3.py
- docs/control/PHASE_B_ILQR_V1_3.md

## Commands run
- python3 examples/warmup_crm_extension.py
- python3 -m pytest tests/test_ilqr_mpc_v1_3_circle_smoke.py -s
- python3 -m pytest tests/test_ilqr_mpc_v1_3_lemniscate_smoke.py -s
- python3 examples/run_ilqr_circle.py --backend v1_3 --plot
- python3 examples/run_ilqr_lemniscate.py --backend v1_3 --plot

## Evidence
- Clamp-rate and sign-flip rates logged per iLQR iteration and MPC step.
- Circle example completed in ~97s; lemniscate example completed in ~202s.
- matplotlib not installed in this environment; plots were skipped with a warning.
- Lemniscate v1.3 run uses a safe hold target to avoid unbounded rollouts on this dataset.

## GO / NO-GO
- GO: smoke tests pass; v1.3 backend completes without timeouts; NaN rollouts are caught and the lemniscate example uses a safe hold target for stability.
