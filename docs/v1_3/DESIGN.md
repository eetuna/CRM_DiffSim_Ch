# v1.3 Implicit Backward Design (Stub)

## Forward equation structure

The CRM dynamics step solves a residual equation of the form:

- G(η, x, u, Li) = 0

where η is the solver decision vector (the nonlinear solver state used internally to satisfy dynamics constraints). The forward output is:

- y = h(η, x, u, Li) = x_{t+1}

The forward step already computes η* and returns x_{t+1} via `crm_step_forward`.

## Implicit differentiation (VJP)

Given a vector v and y = h(η, x, u), we compute vector-Jacobian products without differentiating through solver iterations:

1) Solve the adjoint system:

- (∂G/∂η)^T λ = (∂y/∂η)^T v

2) Then compute:

- v^T dy/du = (∂y/∂u)^T v - (∂G/∂u)^T λ
- v^T dy/dx = (∂y/∂x)^T v - (∂G/∂x)^T λ

No pseudo-inverse, SVD, or damping is allowed in the adjoint solve.

## Definition of η

η is the decision vector used by the CRM nonlinear solver for the dynamics step. It matches the internal solver state used to satisfy the residual equation G(η, x, u) = 0. The exact layout must match the CRM step solver and is cached by `crm_step_forward`.

## Forward cache requirements

Forward must cache (for deterministic backward recomputation):

- η* (solver decision vector)
- solver_exit flag and residual_norm
- any intermediate blocks used to assemble ∂G/∂η, ∂G/∂x, ∂G/∂u, ∂y/∂η, ∂y/∂x, ∂y/∂u
- u0, tau, or other solver outputs needed to reconstruct the residual Jacobians

Cache must be sufficient to recompute the adjoint without re-running solver iterations.

## Linear system dimensions

Let:

- η ∈ R^{n_eta}
- x ∈ R^{n_x}
- u ∈ R^{n_u}
- y ∈ R^{n_y}

Then:

- ∂G/∂η ∈ R^{n_eta × n_eta}
- ∂G/∂x ∈ R^{n_eta × n_x}
- ∂G/∂u ∈ R^{n_eta × n_u}
- ∂y/∂η ∈ R^{n_y × n_eta}
- ∂y/∂x ∈ R^{n_y × n_x}
- ∂y/∂u ∈ R^{n_y × n_u}

Adjoint solve uses transpose conventions as written above.

## Jacobian layout conversion

MINPACK Jacobians are column-major. A single canonical helper converts them to row-major before use:

- `crm_dyn_jacobian_layout_to_row_major(J_col_major)`

All downstream operations consume row-major tensors.

## Non-convergence policy

If forward residual_norm exceeds threshold or solver_exit != 0:

- backward raises an error and reports non-convergence.

No backward is attempted for non-converged forward steps.

## Validation notes (PR-2+)

Full-step finite differences on `crm_step_forward` can be unreliable due to solver path discontinuities.
Correctness is validated via:

- Local FD on blocks: J_eta/J_x/J_u from G(eta, x, u) and H_eta/H_x/H_u from h(eta, x, u).
- Adjoint consistency: check that J_eta^T lambda = H_eta^T v and that the explicit block VJP matches `crm_dyn_vjp_v1_3`.

Current test thresholds:

- Adjoint residual: ||J_eta^T lambda - H_eta^T v|| / (||H_eta^T v|| + 1e-12) < 1e-6
- Block FD comparisons: relative error < 1e-2
- H-block VJP comparisons: max abs diff < 1e-6

## Performance expectations (PR-3)

Implicit VJP avoids autograd Jacobian construction. Expected performance:

- ≥50x speedup versus `torch.autograd.functional.jacobian` on stable inputs.
- Deterministic runtime once the extension is warmed up.

Use v1.3 for iLQR/MPC linearizations whenever you need dense A,B and can remain in the safe regime.
Operator form (A·dx, B·du) is the fastest path and is the preferred API for MPC loops.
Jacobian-free MPC (e.g., CEM/MPPI) remains the fallback when linearization is too expensive or unstable.
