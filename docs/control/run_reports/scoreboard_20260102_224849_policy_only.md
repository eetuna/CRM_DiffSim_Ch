# Run Report (2026-01-02 22:48:50)

- git_commit: 31ca475
- command: examples/scoreboard_eval.py --dataset data/dyn_fk_ramp_circle1_hold1.npz --start-idx 0 --window-len 20 --max-wall-sec 60 --model output_data/models/phase_c2_policy.pt --modes replay_only,policy_only
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
- mean_error: 16.7363395306393
- max_error: 39.29581301916437
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/scoreboard_20260102_224849_policy_only.npz
- plot_prefix: none

- mode: policy_only
- kept_steps: 14
- attempted_steps: 15
- unbounded_count: 0
- solver_exit_counts: {"0": 14, "3": 1}
- runtime_sec: 0.12
- timeout_step: None
