# Phase A iLQR / MPC (CRM differentiable dynamics)

This phase adds a minimal iLQR solver and a receding-horizon MPC wrapper on top of `crm_diffsims.dynamics.crm_step`.

## State and control

- State `x_t`: 39D CRM step state, identical to the `crm_step` input.
  - Tip position used for cost: `x_t[24:27]` (the `xf_pre` tip position block used in CRM dynamics).
- Control `u_t`: actuation currents with shape `[1, 3]` (three coil currents).
- Dynamics: `x_{t+1} = crm_step(x_t, u_t, Li, cfg)` using the existing differentiable step.

## Cost function

Running cost at each step:

- Tip tracking: `w_tip * ||p_tip - p_target||^2`
- Current magnitude: `w_u * ||u_t||^2`
- Current rate: `w_du * ||u_t - u_{t-1}||^2`

Terminal cost:

- Stronger tip tracking penalty: `w_tip_terminal * ||p_tip(T) - p_target(T)||^2`

The rate term is included both in the total cost and in the iLQR backward pass as a local quadratic approximation (including contributions from the adjacent timestep).

## Rate limits and sign-flip prevention

We enforce smoothness during the iLQR forward pass:

- `|u_t - u_{t-1}| <= max_du` clamps current changes each step.
- Large sign flips are limited: if `u_t` crosses zero in one step, the change is capped at `max_sign_flip`.

This prevents solver instability and unrealistic current jumps while keeping the rate penalty in the objective.

## Clamp hit rate diagnostics

Each iLQR solve reports clamp hit rates:

- `clamp_rate`: fraction of timesteps where `|u_t - u_{t-1}|` was clipped to `max_du`.
- `sign_flip_rate`: fraction of timesteps where a sign flip cap was applied.

Frequent clipping indicates a constrained regime where the target or weights demand faster actuation than allowed. If clamp rates are near 100%, reduce target aggressiveness, loosen `max_du`, or increase `w_du` to favor smoother control.

## Tip index diagnostic

`tip_index_sanity_check` logs `x_next[24:27]` for a small control perturbation and verifies it stays finite and changes smoothly. There is no independent tip getter for the step state in the codebase, so this is a lightweight diagnostic to catch obvious indexing mistakes rather than a formal validation.

## Expected operating regime

- Short horizons (20-40 steps) with warm starts between MPC iterations.
- Moderate current magnitudes; aggressive targets may require tighter `max_du` and higher `w_du`.
- Uses CPU-based autograd Jacobians and is intended for validation/debug runs, not high-throughput control.

## Example usage

- `examples/run_ilqr_circle.py`
- `examples/run_ilqr_lemniscate.py`

Each script loads a target trajectory, runs MPC+iLQR, and plots the tip path against the target along with the current profiles.
