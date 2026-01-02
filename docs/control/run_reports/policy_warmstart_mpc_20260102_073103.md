# Run Report (2026-01-02 07:31:03)

- git_commit: 31ca475
- command: examples/eval_policy_warmstart_mpc.py --preset acceptance_fast --max-wall-sec 30 --heartbeat-every 1
- dataset: dyn_fk_ramp_circle1_hold1.npz
- start_idx: 0
- end_idx: 80
- N_total: 0
- Li_mm: 0.0
- dt: 0.0
- umax: 0.0
- d_umax: 0.0
- init_from_ramp: False
- ramp_len: 0
- max_wall_hit: True
- mean_error: nan
- max_error: nan
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/policy_warmstart_20260102_073103.npz
- plot_prefix: none

- timeout_reason: max_wall_sec
- timeout_mode: global
- timeout_step: None
- effective_config: {'preset': 'acceptance_fast', 'start_idx': 0, 'end_idx': 80, 'window_len': 0, 'horizon': 5, 'mpc_steps': 5, 'max_iter': 1, 'max_wall_sec': 60.0, 'preroll_steps': 10, 'init_mode': 'dataset_x0', 'plot': False, 'seed': 0, 'heartbeat_steps': 1}
- eval_go: False
