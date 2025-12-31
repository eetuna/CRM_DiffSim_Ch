# v1.3 (Optional Future Work): True C++ Implicit Backward for Dynamics

This document explains **what v1.3 is**, **why it exists**, **what we have today**, and **why v1.3 is an upgrade—not unfinished work**. It is written to be rigorous but readable and can be used as internal documentation or adapted for a paper appendix.

---

## 1. What we have today (v1.2)

We already have a working, validated differentiable simulator for free-space CRM dynamics.

### Forward dynamics (final and correct)
The forward step is:
\[
x_{t+1} = H(\eta^\*, x_t, u_t)
\]
where the solver finds \(\eta^\*\) by solving:
\[
G(\eta, x_t, u_t) = 0
\]

- This uses the **real legacy physics and solver** (TrustRegionDogleg_dyn / DYNNLEquation).
- The physics is **not approximated**.
- This forward will **not change** in v1.3.

---

## 2. How backward/gradients work in v1.2

### Python finite-difference (FD) backward bridge
In v1.2, gradients are computed using a **Python FD bridge**:

1. Given upstream gradient  
   \[
   g = \frac{\partial L}{\partial x_{t+1}}
   \]
2. The code numerically estimates:
   - \(\left(\frac{\partial x_{t+1}}{\partial \eta}\right)^T g\)
   - \(\left(\frac{\partial x_{t+1}}{\partial u_t}\right)^T g\)
   - \(\left(\frac{\partial x_{t+1}}{\partial x_t}\right)^T g\)
3. This is done by:
   - small perturbations,
   - repeated calls to `crm_step_forward`,
   - finite differencing.

### Important clarification
This is **not “fake gradients.”**

They are:
- mathematically valid,
- deterministic,
- validated by:
  - gradcheck,
  - FD reference comparisons,
  - rollout tests,
  - stress replay tests,
  - workspace branching diagnostics.

v1.2 is therefore **correct and usable** for MPC, iLQR, and RL.

---

## 3. Why the Python FD bridge exists

The legacy dynamics code:
- did **not** expose clean analytic Jacobians for the dynamics residual \(G\),
- relied internally on FD-style Jacobians (MINPACK conventions),
- was written for simulation, not learning.

Jumping directly to a full analytic implicit backward would have required:
- deriving and maintaining many coupled derivatives,
- rewriting solver internals,
- debugging everything at once.

That would have been **high risk and hard to audit**.

The FD bridge allowed us to:
- lock down exact forward semantics,
- define state layout and solver variables precisely,
- enforce determinism and strict failure handling,
- build a comprehensive test and audit harness,
- stress-test initialization and stability behavior.

This was the **correct engineering choice**.

---

## 4. Is the v1.2 backward “fake”?

No. A precise description is:

> **v1.2 has a correct but numerically approximated backward.**

Comparison:

| Property | v1.2 (FD bridge) | v1.3 (implicit backward) |
|--------|------------------|--------------------------|
| Correct forward physics | ✅ | ✅ |
| Deterministic | ✅ | ✅ |
| Valid gradients | ✅ | ✅ |
| Differentiates through solver iterations | ❌ | ❌ |
| Analytic adjoint solve | ❌ | ✅ |
| Large-batch scalability | ⚠️ | ✅ |
| Numerical noise | ⚠️ | ✅ |
| Performance | ⚠️ | ✅ |

---

## 5. What v1.3 actually changes

### v1.3 does **not** change:
- the forward dynamics,
- the solver,
- the physics,
- the public API.

### v1.3 **does** change:
- how gradients are computed,
- where they are computed (C++ instead of Python),
- performance and numerical smoothness.

---

## 6. The correct implicit backward (v1.3 mathematics)

Define:
- \(A = \partial G / \partial \eta\)
- \(B = \partial G / \partial u_t\)
- \(C = \partial G / \partial x_t\)

and:
- \(H_\eta = \partial H / \partial \eta\)
- \(H_u = \partial H / \partial u_t\)
- \(H_x = \partial H / \partial x_t\)

Given upstream gradient \(g = \partial L / \partial x_{t+1}\):

1. **Adjoint solve**:
\[
A^T \lambda = H_\eta^T g
\]

2. **Gradient computation**:
\[
\frac{\partial L}{\partial u_t}
=
H_u^T g - B^T \lambda
\]

\[
\frac{\partial L}{\partial x_t}
=
H_x^T g - C^T \lambda
\]

This is the **same implicit-differentiation structure used in v1.0 (equilibrium)**, now applied to dynamics.

---

## 7. What code v1.3 would add

- A C++ implementation of `crm_step_backward_impl(...)`.
- Assembly of:
  - `A_step`, `B_step`, `C_step`,
  - `H_eta^T g`, `H_u^T g`, `H_x^T g`.
- One **small (6×6) linear solve per step**.
- Direct return of gradients to PyTorch.

### What disappears
- Python FD loops,
- repeated forward calls during backward,
- FD step-size tuning,
- FD-induced gradient noise.

---

## 8. Why v1.3 is better (when you need it)

### Performance
- FD backward: many forward calls per gradient.
- Implicit backward: one small linear solve.
- Critical for:
  - long horizons,
  - large batches,
  - RL training loops,
  - repeated MPC solves.

### Numerical quality
- FD gradients are step-size sensitive and noisy near singularities.
- Implicit gradients are smoother and more stable.

### Scientific clarity
With v1.3 you can state:
> “Dynamics gradients are computed via analytic implicit differentiation without differentiating through solver iterations.”

---

## 9. Can v1.1, v1.2, and v1.3 coexist?

Yes.

- v1.1/v1.2 define the **semantics and API** of `crm_step`.
- v1.3 is an **internal upgrade**: same forward, same API, better backward.
- FD reference paths can remain behind debug flags for validation.

---

## 10. Why v1.3 is new scope, not unfinished work

v1.2 already provides:
- correct forward physics,
- validated gradients,
- deterministic behavior,
- strict failure policy,
- stress-tested stability and initialization,
- reproducible audit trail.

v1.3 is a **performance and rigor upgrade**, not a correctness fix.

---

## 11. When v1.3 is worth doing

Pursue v1.3 if you need:
- large-batch RL,
- long-horizon backpropagation,
- faster training,
- publication-level claims of fully analytic implicit dynamics gradients.

If not, **v1.2 is already a strong and defensible endpoint**.

---

## 12. One-sentence summary

**v1.2 already works and is validated; v1.3 replaces the Python finite-difference backward with a true C++ implicit adjoint backward, improving speed and gradient quality while leaving the forward dynamics unchanged.**
