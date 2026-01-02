# Run Report (2026-01-02 22:52:43)

- git_commit: 31ca475
- command: examples/scoreboard_eval.py --dataset data/dyn_fk_ramp_circle1_hold1.npz --start-idx 0 --window-len 50 --modes replay_only,policy_only,policy_plus_mpc --max-wall-sec 60
- dataset: dyn_fk_ramp_circle1_hold1.npz
- start_idx: 0
- end_idx: 50
- N_total: 320
- Li_mm: 94.3
- dt: 0.05
- umax: 0.1
- d_umax: 0.01
- init_from_ramp: False
- ramp_len: 0
- max_wall_hit: True
- mean_error: 38.23269043363562
- max_error: 39.29581301916437
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/scoreboard_20260102_225133_policy_plus_mpc.npz
- plot_prefix: none

- mode: policy_plus_mpc
- kept_steps: 2
- attempted_steps: 2
- unbounded_count: 0
- solver_exit_counts: {"0": 2}
- runtime_sec: 69.44
- timeout_step: 2
