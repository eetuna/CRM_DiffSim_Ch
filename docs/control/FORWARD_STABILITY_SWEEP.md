# Forward Stability Sweep

This sweep runs forward-only CRM rollouts (no Jacobians) to identify when NaNs or "Unbounded" warnings appear as current magnitude and rate limits increase.

## Usage

- Basic:
  - `python3 examples/sweep_forward_stability.py --umax_list 0.05,0.1,0.2 --d_umax_list 0.01,0.02 --T 50`
- With dataset currents:
  - `python3 examples/sweep_forward_stability.py --dataset data/dyn_fk_ramp_circle1_hold1.npz --umax_list 0.05,0.1 --d_umax_list 0.01 --T 50`

## Outputs

- CSV: `output_data/forward_stability_sweep.csv`
- Failures: `output_data/forward_failures/flight_*.npz` and `.json`

## Interpretation

- `nan_detected`: forward dynamics produced NaN/Inf.
- `unbounded_detected`: "Coil integration Unbounded!!" appears in captured output (note: C++ output may not always be captured).
- `mean_step_time` and `max_step_time`: forward performance for the chosen (umax, d_umax) pair.
- `fail_index`: step index where NaN/Inf first appears (or -1).

If failures only appear at larger `umax` or `d_umax`, the instability is likely due to control magnitude/rate rather than baseline dynamics.
