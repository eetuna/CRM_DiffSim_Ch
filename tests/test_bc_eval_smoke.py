import os
import subprocess
import sys


def test_bc_eval_smoke():
    dataset_npz = "output_data/test_bc_eval_dataset.npz"
    dataset_report = "docs/control/run_reports/test_bc_eval_dataset.md"
    model_path = "output_data/models/test_bc_eval_policy.pt"
    train_npz = "output_data/test_bc_eval_train_metrics.npz"
    train_report = "docs/control/run_reports/test_bc_eval_train_report.md"

    for path in [dataset_npz, dataset_report, model_path, train_npz, train_report]:
        if os.path.exists(path):
            os.remove(path)

    gen_cmd = [
        sys.executable,
        "examples/generate_bc_dataset.py",
        "--controller",
        "replay_only_full",
        "--episodes",
        "1",
        "--episode-len",
        "8",
        "--start-idx",
        "0",
        "--end-idx",
        "16",
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
        "1",
        "--batch-size",
        "4",
        "--hidden",
        "16",
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

    eval_prefix = "test_bc_eval"
    report_prefix = "test_bc_eval_report"
    eval_cmd = [
        sys.executable,
        "examples/eval_bc_policy_rollout.py",
        "--model",
        model_path,
        "--start-idx",
        "0",
        "--subset-len",
        "6",
        "--output-prefix",
        eval_prefix,
        "--report-prefix",
        report_prefix,
    ]
    subprocess.run(eval_cmd, check=True)

    circle_npz = f"output_data/{eval_prefix}_circle.npz"
    lemniscate_npz = f"output_data/{eval_prefix}_lemniscate.npz"
    circle_report = f"docs/control/run_reports/{report_prefix}_circle.md"
    lemniscate_report = f"docs/control/run_reports/{report_prefix}_lemniscate.md"

    assert os.path.exists(circle_npz)
    assert os.path.exists(lemniscate_npz)
    assert os.path.exists(circle_report)
    assert os.path.exists(lemniscate_report)
