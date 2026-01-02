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
