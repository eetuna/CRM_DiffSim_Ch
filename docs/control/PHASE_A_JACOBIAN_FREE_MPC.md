# Phase A Jacobian-Free MPC (CEM)

Autograd Jacobians for the CRM step are prohibitively expensive for full-fidelity MPC. This baseline uses forward-only rollouts with a Cross-Entropy Method (CEM) shooting optimizer.

## Why Jacobian-free

- Forward step is ~milliseconds; Jacobians are orders of magnitude slower.
- CEM avoids Jacobians entirely by sampling control sequences and scoring rollouts.

## Algorithm sketch

- Sample N candidate control sequences around a mean (with diagonal std).
- Clamp each candidate to safe bounds (|u|, |Δu|, optional sign-flip cap).
- Roll out dynamics and compute cost.
- Keep top-K elites, refit mean/std.
- Apply the first control, shift the mean for warm start.

## Safe bounds

All candidates and applied controls are hard-clamped to `safe_bounds.json` (derived from sweep) or defaults if missing. This keeps the controller in a regime that avoids unbounded coil warnings and NaNs.

## Plotting and Inspection

Use the `--plot` flag on the CEM examples to save high-resolution PNG and PDF plots under `output_data/`. Plotting is optional and requires the `plots` extra (`matplotlib`). If matplotlib is unavailable, the examples print a warning and continue without plotting.
