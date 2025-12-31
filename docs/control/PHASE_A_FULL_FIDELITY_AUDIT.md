# Phase A Full-Fidelity Audit

## Reproduction commands

- `python3 examples/warmup_crm_extension.py`
- `python3 -m pytest tests/test_forward_rollout_stability.py -s`
- `python3 -m pytest tests/test_linearization_cost.py -s`
- `python3 examples/run_ilqr_circle.py --profile-step --max-wall-sec 300 --save-u-seq`
- `python3 examples/run_ilqr_lemniscate.py --profile-step --max-wall-sec 300 --save-u-seq`

## Measurements

- Forward rollout stability (T=20 ramp, no optimization):
  - mean_step_time ~0.0131s
  - no NaN/Inf, no "Unbounded" warnings
- Full-fidelity linearization cost:
  - step_time ~0.0107s
  - linearization_time ~134.38s (1 step)
  - ratio ~12508x (linearization dominates)
- Partial full-fidelity MPC runs (max-wall 300s):
  - circle profile mean_step ~0.2149s for the applied step
  - lemniscate profile mean_step ~0.3592s for the applied step
  - both exit due to wall time before completing MPC horizon

## Failure characterization

- **Forward instability (A)**:
  - Not observed in the forward-only stability test with small smooth currents.
  - "Coil integration Unbounded!!" did not appear in this test.
- **Linearization cost (B)**:
  - Dominant. Autograd Jacobian computation is ~10^4x the forward step time.
  - Full-fidelity iLQR spends most time inside `_linearize`.
- **MPC configuration (C)**:
  - Long horizons (20-40) and multiple iLQR iterations multiply the Jacobian cost.
  - Without a trust-region cut or very short horizons, the wall time exceeds 300s quickly.
- **Unbounded warnings + NaNs**:
  - Unbounded warnings previously appeared in lemniscate full runs; not reproduced in the forward-only ramp.
  - Likely triggered by large control magnitudes during optimization rather than base dynamics.

## Conclusion

Root cause: **Linearization cost (B)** is the primary reason full-fidelity runs are extremely slow. Forward dynamics are stable under small smooth currents (A not dominant). MPC configuration (C) exacerbates the cost by multiplying Jacobian evaluations across the horizon and iterations. The "Unbounded" warnings are likely triggered by large currents during optimization but do not appear in the forward-only stability test.

## Recommended guardrails and parameter adjustments

- Reduce horizon for full-fidelity debug runs: T <= 10.
- Reduce iLQR iterations for full-fidelity runs: max_iter <= 3.
- Increase `w_du`, reduce target aggressiveness to prevent large u.
- Clamp currents via `max_u` to keep coil integration stable.
- Use `--max-wall-sec` and `--profile-step` in examples to capture partial traces and diagnose stalls.
- Consider caching/parallelizing Jacobian computations if/when permitted (outside Phase A constraints).
