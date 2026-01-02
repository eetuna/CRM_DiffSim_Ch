# Run Report Contract

Every debug/run must generate:
1) A run report under `docs/control/run_reports/`
2) An NPZ artifact under `output_data/`
3) PNG+PDF plots when `--plot` is set

The run report must include:
- Command line used
- Git commit hash (if available)
- Dataset name and slice (`start_idx`, `end_idx`, `N_total`)
- Li and `dt`
- Bounds used (`umax`, `d_umax`)
- Ramp settings (`init-from-ramp`, `ramp_len`)
- Mean/max error
- Max-wall hit status
- Artifact path and plot prefix
- Optional diagnostics (when available):
  - `requested_start_idx`, `requested_end_idx`
  - `ext_load_sec`, `dataset_load_sec`, `model_load_sec`, `rollout_elapsed_sec`, `total_elapsed_sec`
  - `mpc_only_exception_type`, `mpc_only_exception_msg`, `mpc_only_exception_tb`
