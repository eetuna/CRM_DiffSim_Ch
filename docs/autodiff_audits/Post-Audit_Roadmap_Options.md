# Post-Audit Roadmap Options (A–F)

This document records all valid next steps after completing:
- v1.0 / v1.1 / v1.2 implementation
- full gradient validation
- stress testing and initialization audit
- audit close-out and branch freeze

All options below are **independent** and **optional**. None represent unfinished work.

---

## Option A — USE THE SIMULATOR (Primary / Default)
**Status:** Ready now  
**Risk:** Low  
**Value:** High

### Description
Begin *using* the differentiable CRM simulator as a tool for:
- MPC / iLQR
- trajectory optimization
- model-based RL
- residual learning
- offline or online evaluation

### Why this is valid
The simulator already provides:
- correct forward physics
- validated gradients
- deterministic behavior
- stress-tested stability
- understood initialization constraints
- reproducible audit trail

No further simulator engineering is required to do real work.

### Typical tasks
- Define cost functions (tip error, current smoothness, constraints)
- Implement an MPC or iLQR loop
- Run rollouts on workspace trajectories
- Analyze performance and convergence

---

## Option B — PAPER / METHODS / DOCUMENTATION
**Status:** Ready now  
**Risk:** Low  
**Value:** Very high

### Description
Convert the completed work into:
- a Methods section
- a Supplementary audit appendix
- a standalone systems paper

### Why now is ideal
- All design decisions are fresh
- Stress tests and initialization analysis are rare and publishable
- Full reproducibility already exists

v1.2 is more than sufficient for publication.

---

## Option C — v1.3: TRUE C++ IMPLICIT BACKWARD (Future Research)
**Status:** New scope  
**Risk:** Medium  
**Value:** High (conditional)

### Description
Replace the Python finite-difference backward with a **true C++ implicit adjoint backward**:
- same forward dynamics
- same API
- faster and cleaner gradients

### What this is NOT
- Not new physics
- Not re-deriving Cosserat theory
- Not unfinished work

### When this is worth doing
Only if real usage (Option A) shows:
- FD backward is a bottleneck
- large-batch RL is needed
- long-horizon gradients are required

---

## Option D — HARDWARE / MRI / EXPERIMENTAL VALIDATION
**Status:** Depends on setup  
**Risk:** Medium  
**Value:** Very high

### Description
Use the simulator as a digital twin to:
- replay experimental current traces
- compare predicted vs measured tip motion
- tune parameters
- validate under MRI constraints

Strongly aligned with MRI-guided catheter research.

---

## Option E — RESIDUAL / LEARNED MODEL ON TOP OF CRM
**Status:** Ready  
**Risk:** Medium  
**Value:** High

### Description
Keep CRM as the physics backbone and learn:
- residual forces
- unmodeled damping
- hysteresis
- magnetic field imperfections

A common and effective hybrid modeling strategy.

---

## Option F — MULTI-ACTUATOR / CONTACT / CONSTRAINTS
**Status:** Large new scope  
**Risk:** High  
**Value:** Very high

### Description
Extend beyond free-space single-actuator dynamics:
- multiple actuators
- vessel wall contact
- friction
- contact-implicit dynamics

This is **new research**, not continuation.

---

## Recommended Order
**Option A → Option B → (only if needed) Option C**

Use the simulator first, document it properly, then optimize if justified.

---

## Mental Reset
All milestones are closed.
Everything from here is **choice**, not obligation.
