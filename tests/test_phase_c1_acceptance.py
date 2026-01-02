import json
import os
import subprocess
import sys

import numpy as np


def _read_drop_reasons(path):
    with np.load(path, allow_pickle=True) as data:
        reasons = data.get("drop_reasons_counts", "{}")
        if isinstance(reasons, np.ndarray):
            reasons = reasons.item()
        if isinstance(reasons, bytes):
            reasons = reasons.decode("ascii")
        try:
            reasons = json.loads(reasons)
        except Exception:
            reasons = {}
    return reasons


def test_phase_c1_acceptance():
    dataset_npz = "output_data/phase_c1_dataset.npz"
    dataset_report = "docs/control/run_reports/phase_c1_dataset.md"
    model_path = "output_data/models/phase_c1_policy.pt"
    train_npz = "output_data/phase_c1_train_metrics.npz"
    train_report = "docs/control/run_reports/phase_c1_train_report.md"
    eval_prefix = "phase_c1_eval"
    eval_report_prefix = "phase_c1_eval_report"

    for path in [dataset_npz, dataset_report, model_path, train_npz, train_report]:
        if os.path.exists(path):
            os.remove(path)

    gen_cmd = [
        sys.executable,
        "examples/generate_bc_dataset.py",
        "--preset",
        "mixed_easy",
        "--episodes",
        "2",
        "--episode-len",
        "30",
        "--filter-policy",
        "drop_unbounded",
        "--seed",
        "0",
        "--output-npz",
        dataset_npz,
        "--report-path",
        dataset_report,
    ]
    subprocess.run(gen_cmd, check=True)

    train_cmd = [
        sys.executable,
        "examples/train_bc_policy.py",
        "--dataset",
        dataset_npz,
        "--epochs",
        "3",
        "--batch-size",
        "8",
        "--hidden",
        "32",
        "--include-tip",
        "--seed",
        "0",
        "--output-model",
        model_path,
        "--output-npz",
        train_npz,
        "--report-path",
        train_report,
    ]
    subprocess.run(train_cmd, check=True)

    eval_cmd = [
        sys.executable,
        "examples/eval_bc_policy_rollout.py",
        "--model",
        model_path,
        "--start-idx",
        "0",
        "--subset-len",
        "20",
        "--output-prefix",
        eval_prefix,
        "--report-prefix",
        eval_report_prefix,
        "--bc-dataset",
        dataset_npz,
    ]
    subprocess.run(eval_cmd, check=True)

    with np.load(dataset_npz, allow_pickle=True) as data:
        total = int(data.get("n_steps_total", 0))
        kept = int(data.get("n_steps_kept", 0))
        dropped = int(data.get("n_steps_dropped", 0))

    frac_dropped = dropped / max(total, 1)
    print(f"phase_c1: dropped_fraction={frac_dropped:.3f} ({dropped}/{total})")
    reasons = _read_drop_reasons(dataset_npz)
    print(f"phase_c1: drop_reasons={reasons}")
    assert frac_dropped <= 0.30

    circle_npz = f"output_data/{eval_prefix}_circle.npz"
    lem_npz = f"output_data/{eval_prefix}_lemniscate.npz"
    assert os.path.exists(circle_npz)
    assert os.path.exists(lem_npz)

    with np.load(circle_npz) as data:
        err_pol = data["error_policy"]
        assert np.isfinite(err_pol).all()
        mean_err = float(np.mean(err_pol)) if err_pol.size else 0.0
        assert mean_err < 50.0

    with np.load(lem_npz) as data:
        err_pol = data["error_policy"]
        assert np.isfinite(err_pol).all()
