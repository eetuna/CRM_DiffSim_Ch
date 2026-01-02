import os
import subprocess
import sys

import numpy as np


def test_bc_eval_keeps_steps_real_subset():
    dataset_npz = "output_data/test_eval_subset_dataset.npz"
    dataset_report = "docs/control/run_reports/test_eval_subset_dataset.md"
    model_path = "output_data/models/test_eval_subset_policy.pt"
    train_npz = "output_data/test_eval_subset_train.npz"
    train_report = "docs/control/run_reports/test_eval_subset_train.md"
    eval_prefix = "test_eval_subset"
    report_prefix = "test_eval_subset_report"

    for path in [dataset_npz, dataset_report, model_path, train_npz, train_report]:
        if os.path.exists(path):
            os.remove(path)

    gen_cmd = [
        sys.executable,
        "examples/generate_bc_dataset.py",
        "--preset",
        "circle_easy",
        "--episodes",
        "1",
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
        "2",
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
        "--init-mode",
        "dataset_x0",
        "--start-idx",
        "0",
        "--subset-len",
        "20",
        "--output-prefix",
        eval_prefix,
        "--report-prefix",
        report_prefix,
    ]
    subprocess.run(eval_cmd, check=True)

    circle_npz = f"output_data/{eval_prefix}_circle.npz"
    assert os.path.exists(circle_npz)

    with np.load(circle_npz, allow_pickle=True) as data:
        keep = data["policy_keep_mask"].astype(bool)
        tip = data["tip_policy"]
        kept_steps = int(keep.sum())
        assert kept_steps >= 10
        if kept_steps:
            assert np.isfinite(tip[keep]).all()
