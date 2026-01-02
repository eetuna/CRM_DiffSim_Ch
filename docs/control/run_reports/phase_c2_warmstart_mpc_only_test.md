# Run Report (2026-01-02 21:45:52)

- git_commit: 31ca475
- command: examples/eval_policy_warmstart_mpc.py --model output_data/models/phase_c2_policy.pt --dataset data/dyn_fk_ramp_circle1_hold1.npz --start-idx 0 --window-len 50 --init-mode dataset_x0 --preroll-steps 5 --horizon 5 --mpc-steps 3 --max-wall-sec 120 --seed 0 --output-prefix phase_c2_warmstart_mpc_only_test --report-path docs/control/run_reports/phase_c2_warmstart_mpc_only_test.md
- dataset: dyn_fk_ramp_circle1_hold1.npz
- start_idx: 5
- end_idx: 55
- N_total: 320
- Li_mm: 94.3
- dt: 0.05
- umax: 0.1
- d_umax: 0.01
- init_from_ramp: False
- ramp_len: 0
- max_wall_hit: False
- mean_error: 7.195470516469247
- max_error: 15.34742372299618
- baseline_mean_error: n/a
- baseline_max_error: n/a
- artifact_npz: output_data/phase_c2_warmstart_mpc_only_test.npz
- plot_prefix: none

- requested_start_idx: 0
- requested_end_idx: 50
- init_mode: dataset_x0
- prev_u_source: preroll_u_last
- base_start_idx: 0
- preroll_steps: 5
- preroll_start_idx: 0
- preroll_end_idx: 5
- max_wall_sec: 120.0
- mpc_steps: 3
- ext_load_sec: 1.053415060043335
- dataset_load_sec: 0.005134105682373047
- model_load_sec: 0.00781702995300293
- rollout_elapsed_sec: 57.51134729385376
- total_elapsed_sec: 58.757001638412476
- effective_config: {'preset': 'none', 'start_idx': 5, 'end_idx': 55, 'window_len': 50, 'horizon': 5, 'mpc_steps': 3, 'max_iter': 1, 'max_wall_sec': 120.0, 'preroll_steps': 5, 'init_mode': 'dataset_x0', 'plot': False, 'seed': 0, 'heartbeat_steps': 10}
- policy_only_mean: 7.508401838944771
- policy_only_max: 15.34742372299618
- mpc_only_mean: 10.005735862123784
- mpc_only_max: 15.34742372299618
- policy_plus_mpc_mean: 7.195470516469247
- policy_plus_mpc_max: 15.34742372299618
- policy_only_kept: 5
- mpc_only_kept: 3
- policy_plus_mpc_kept: 9
- policy_only_runtime_sec: 0.08996343612670898
- mpc_only_runtime_sec: 9.490966081619263
- policy_plus_mpc_runtime_sec: 47.963146448135376
- policy_only_timeout: False
- mpc_only_timeout: False
- policy_plus_mpc_timeout: False
- policy_only_timeout_step: None
- mpc_only_timeout_step: None
- policy_plus_mpc_timeout_step: None
- policy_only_timeout_reason: none
- mpc_only_timeout_reason: none
- policy_plus_mpc_timeout_reason: none
- mpc_only_exception_type: IndexError
- mpc_only_exception_msg: index 1 is out of bounds for dimension 0 with size 1
- mpc_only_exception_tb: Traceback (most recent call last):
  File "/workspaces/CRM_DiffSim_Ch/examples/eval_policy_warmstart_mpc.py", line 454, in _mpc_loop
    result = ilqr_solve(
  File "/workspaces/CRM_DiffSim_Ch/python/crm_diffsims/control/ilqr.py", line 282, in ilqr_solve
    cost = total_cost(tip_seq, u_seq, targets_local, cfg)
  File "/workspaces/CRM_DiffSim_Ch/python/crm_diffsims/control/ilqr.py", line 123, in total_cost
    tip = tip_seq[t + 1]
IndexError: index 1 is out of bounds for dimension 0 with size 1
- failure_reason_mpc_only: 
- first_fail_step_mpc_only: n/a
- first_fail_solver_exit_mpc_only: n/a
- first_fail_residual_mpc_only: n/a
- first_fail_unbounded_mpc_only: n/a
- eval_go: True
