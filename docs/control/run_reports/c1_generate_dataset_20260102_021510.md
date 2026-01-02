# Run Report (2026-01-02 02:20:10)

- git_commit: 31ca475
- command: examples/generate_bc_dataset.py --preset mixed_easy --episodes 240 --episode-len 400 --wrap-dataset --filter-policy drop_unbounded --seed 0 --output-npz output_data/c1_bc_dataset_20260102_021510.npz --report-path docs/control/run_reports/c1_generate_dataset_20260102_021510.md
- dataset: dyn_fk_ramp_circle1_hold1.npz,dyn_fk_lem1_y40_a10_L94_hold1.npz
- start_idx: 0
- end_idx: -1
- N_total: 640
- Li_mm: 94.3
- dt: 0.05
- umax: 0.1
- d_umax: 0.01
- init_from_ramp: True
- ramp_len: 120
- max_wall_hit: False
- mean_error: nan
- max_error: nan
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/c1_bc_dataset_20260102_021510.npz
- plot_prefix: none

- filter_policy: drop_unbounded
- n_steps_total: 12780
- n_steps_kept: 12660
- n_steps_dropped: 120
- drop_reasons_counts: {"unbounded": 0, "solver_exit": 60, "residual": 120, "nonfinite": 0, "policy_drop": 0, "prefix_fail": 120}
