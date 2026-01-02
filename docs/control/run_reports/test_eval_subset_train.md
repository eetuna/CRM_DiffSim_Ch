# Run Report (2026-01-02 03:12:30)

- git_commit: 31ca475
- command: examples/train_bc_policy.py --dataset output_data/test_eval_subset_dataset.npz --epochs 2 --batch-size 8 --hidden 32 --include-tip --seed 0 --output-model output_data/models/test_eval_subset_policy.pt --output-npz output_data/test_eval_subset_train.npz --report-path docs/control/run_reports/test_eval_subset_train.md
- dataset: dyn_fk_ramp_circle1_hold1.npz
- start_idx: 0
- end_idx: -1
- N_total: 320
- Li_mm: 94.3
- dt: 0.05
- umax: 0.1
- d_umax: 0.01
- init_from_ramp: True
- ramp_len: 120
- max_wall_hit: False
- mean_error: 1.5723658800125122
- max_error: 2.0040769577026367
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/test_eval_subset_train.npz
- plot_prefix: none

- n_steps_total: 30
- n_steps_kept: 30
- n_steps_dropped: 0
- drop_reasons_counts: {'unbounded': 0, 'solver_exit': 0, 'residual': 0, 'nonfinite': 0, 'target_nonfinite': 0, 'policy_drop': 0, 'prefix_fail': 0}
