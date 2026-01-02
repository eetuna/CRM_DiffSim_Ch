# Run Report (2026-01-02 22:51:20)

- git_commit: 31ca475
- command: examples/generate_bc_dataset.py --controller cem_mpc --episodes 1 --episode-len 20 --start-idx 0 --end-idx 20 --horizon 5 --mpc-steps 5 --num-samples 8 --num-elites 2 --cem-iters 1 --seed 0 --output-npz output_data/test_bc_dataset_cem_smoke.npz --report-path docs/control/run_reports/test_bc_dataset_cem_smoke.md
- dataset: dyn_fk_ramp_circle1_hold1.npz
- start_idx: 0
- end_idx: 20
- N_total: 320
- Li_mm: 94.3
- dt: 0.05
- umax: 0.1
- d_umax: 0.01
- init_from_ramp: False
- ramp_len: 0
- max_wall_hit: False
- mean_error: 9.829316281730735
- max_error: 39.29581301916437
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/test_bc_dataset_cem_smoke.npz
- plot_prefix: none

- controller: cem_mpc
- filter_policy: drop_unbounded
- n_steps_total: 20
- n_steps_kept: 20
- n_steps_dropped: 0
- drop_reasons_counts: {"unbounded": 0, "solver_exit": 0, "residual": 0, "nonfinite": 0, "target_nonfinite": 0, "policy_drop": 0, "prefix_fail": 0}
