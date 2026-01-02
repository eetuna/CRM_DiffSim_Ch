# Run Report (2026-01-02 22:22:51)

- git_commit: 31ca475
- command: examples/eval_policy_warmstart_mpc.py --model output_data/models/phase_c2_policy.pt --dataset data/dyn_fk_ramp_circle1_hold1.npz --start-idx 0 --window-len 50 --init-mode dataset_x0 --preroll-steps 5 --horizon 5 --mpc-steps 3 --max-wall-sec 300 --heartbeat-every 1 --seed 0 --output-prefix phase_c2_warmstart_debug --report-path docs/control/run_reports/phase_c2_warmstart_debug.md
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
- artifact_npz: output_data/phase_c2_warmstart_debug.npz
- plot_prefix: none

- requested_start_idx: 0
- requested_end_idx: 50
- init_mode: dataset_x0
- prev_u_source: preroll_u_last
- base_start_idx: 0
- preroll_steps: 5
- preroll_start_idx: 0
- preroll_end_idx: 5
- max_wall_sec: 300.0
- mpc_steps: 3
- ext_load_sec: 1.1550326347351074
- dataset_load_sec: 0.014290332794189453
- model_load_sec: 0.016346454620361328
- rollout_elapsed_sec: 49.76856517791748
- total_elapsed_sec: 51.16141867637634
- effective_config: {'preset': 'none', 'start_idx': 5, 'end_idx': 55, 'window_len': 50, 'horizon': 5, 'mpc_steps': 3, 'max_iter': 1, 'max_wall_sec': 300.0, 'preroll_steps': 5, 'init_mode': 'dataset_x0', 'plot': False, 'seed': 0, 'heartbeat_steps': 1}
- policy_only_mean: 7.508401838944771
- policy_only_max: 15.34742372299618
- mpc_only_mean: 10.005735862123784
- mpc_only_max: 15.34742372299618
- policy_plus_mpc_mean: 7.195470516469247
- policy_plus_mpc_max: 15.34742372299618
- policy_only_kept: 5
- mpc_only_kept: 3
- policy_plus_mpc_kept: 9
- policy_only_runtime_sec: 0.07356834411621094
- mpc_only_runtime_sec: 9.743047952651978
- policy_plus_mpc_runtime_sec: 39.9859983921051
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
