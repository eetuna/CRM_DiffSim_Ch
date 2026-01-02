import os
import subprocess
from datetime import datetime


def _git_hash():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return "unknown"
    return result.stdout.strip() or "unknown"


def write_run_report(
    *,
    report_path,
    command,
    dataset_name,
    start_idx,
    end_idx,
    n_total,
    li_mm,
    dt,
    umax,
    d_umax,
    mean_error,
    max_error,
    artifact_npz,
    plot_prefix,
    init_from_ramp=False,
    ramp_len=0,
    max_wall_hit=False,
    baseline_mean=None,
    baseline_max=None,
):
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    git_hash = _git_hash()
    lines = [
        f"# Run Report ({timestamp})",
        "",
        f"- git_commit: {git_hash}",
        f"- command: {command}",
        f"- dataset: {dataset_name}",
        f"- start_idx: {start_idx}",
        f"- end_idx: {end_idx}",
        f"- N_total: {n_total}",
        f"- Li_mm: {li_mm}",
        f"- dt: {dt}",
        f"- umax: {umax}",
        f"- d_umax: {d_umax}",
        f"- init_from_ramp: {bool(init_from_ramp)}",
        f"- ramp_len: {ramp_len}",
        f"- max_wall_hit: {bool(max_wall_hit)}",
        f"- mean_error: {mean_error}",
        f"- max_error: {max_error}",
        f"- baseline_mean_error: {baseline_mean if baseline_mean is not None else 'n/a'}",
        f"- baseline_max_error: {baseline_max if baseline_max is not None else 'n/a'}",
        f"- artifact_npz: {artifact_npz}",
        f"- plot_prefix: {plot_prefix if plot_prefix is not None else 'none'}",
        "",
    ]
    with open(report_path, "w", encoding="ascii") as f:
        f.write("\n".join(lines))
