# Phase C1 Acceptance

## Dataset hygiene
- source: mixed_easy preset
- steps_total: 60
- steps_kept: 60
- steps_dropped: 0
- drop_reasons_counts: {"unbounded": 0, "solver_exit": 0, "residual": 0, "nonfinite": 0, "target_nonfinite": 0, "policy_drop": 0, "prefix_fail": 0}

## Training
- loss_initial: 8.4171082620
- loss_final: 1.1201585884

## Evaluation
- circle: mean_error=12.8178056485 max_error=39.2958130192
- lemniscate: mean_error=13.3796809366 max_error=39.2958130192
- policy_unbounded_frac: 0.0

## Decision
GO for Phase C2 (MPC-in-the-loop residual / warm-start policy).

## Criteria update
- Eval must keep at least one step (`kept_steps > 0`); otherwise report NO-GO.
