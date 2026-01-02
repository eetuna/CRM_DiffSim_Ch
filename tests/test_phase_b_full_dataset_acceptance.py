import glob
import os
import subprocess
import sys
from datetime import datetime

import numpy as np

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
python_root = os.path.join(repo_root, "python")
if python_root not in sys.path:
    sys.path.insert(0, python_root)

from crm_diffsims.control.run_report import write_run_report


def _run_and_get_npz(cmd, pattern, allow_timeout=False):
    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    before = set(glob.glob(os.path.join(out_dir, pattern)))
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
    if allow_timeout:
        assert result.returncode in (0, 2)
    else:
        assert result.returncode == 0
    after = set(glob.glob(os.path.join(out_dir, pattern)))
    new_files = sorted(after - before)
    assert new_files
    return new_files[-1]


def _mean_max_error(npz_path):
    data = np.load(npz_path, allow_pickle=False)
    tip = data["executed_tip_xyz"]
    target = data["target_xyz"]
    n = min(tip.shape[0], target.shape[0])
    errors = np.linalg.norm(tip[:n] - target[:n], axis=1)
    return float(errors.mean()), float(errors.max()), data


def test_phase_b_full_dataset_acceptance():
    report_dir = os.path.join("docs", "control", "run_reports")
    os.makedirs(report_dir, exist_ok=True)
    report_before = set(glob.glob(os.path.join(report_dir, "*.md")))

    baseline_cmd = [
        sys.executable,
        "examples/run_ilqr_circle.py",
        "--backend",
        "v1_3",
        "--full-dataset",
        "--mode",
        "replay_only_full",
        "--start-idx",
        "0",
        "--end-idx",
        "150",
        "--max-wall-sec",
        "300",
    ]
    baseline_npz = _run_and_get_npz(baseline_cmd, "replay_mode_circle_*.npz")
    baseline_mean, baseline_max, baseline_data = _mean_max_error(baseline_npz)

    ctrl_cmd = [
        sys.executable,
        "examples/run_ilqr_circle.py",
        "--backend",
        "v1_3",
        "--full-dataset",
        "--mode",
        "replay_dataset_mode_full",
        "--start-idx",
        "0",
        "--end-idx",
        "150",
        "--init-from-ramp",
        "--ramp-len",
        "50",
        "--max-iter",
        "1",
        "--max-wall-sec",
        "300",
    ]
    ctrl_npz = _run_and_get_npz(ctrl_cmd, "replay_mode_circle_*.npz", allow_timeout=True)
    ctrl_mean, ctrl_max, ctrl_data = _mean_max_error(ctrl_npz)

    assert np.isfinite(ctrl_data["executed_tip_xyz"]).all()
    assert ctrl_mean <= baseline_mean + 10.0

    report_after = set(glob.glob(os.path.join(report_dir, "*.md")))
    assert report_after - report_before

    comparison_path = os.path.join(
        report_dir, f"phase_b_acceptance_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
    )
    write_run_report(
        report_path=comparison_path,
        command="; ".join([" ".join(baseline_cmd), " ".join(ctrl_cmd)]),
        dataset_name="data/dyn_fk_ramp_circle1_hold1.npz",
        start_idx=0,
        end_idx=150,
        n_total=150,
        li_mm=float(ctrl_data["Li_mm"]),
        dt=float(ctrl_data["dt"]),
        umax=float(ctrl_data["umax_safe"]),
        d_umax=float(ctrl_data["d_umax_safe"]),
        mean_error=ctrl_mean,
        max_error=ctrl_max,
        artifact_npz=ctrl_npz,
        plot_prefix=None,
        init_from_ramp=True,
        ramp_len=50,
        baseline_mean=baseline_mean,
        baseline_max=baseline_max,
    )
    assert os.path.exists(comparison_path)
