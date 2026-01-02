import os
import subprocess
import sys


def test_bc_train_smoke():
    dataset_npz = "output_data/test_bc_train_dataset.npz"
    dataset_report = "docs/control/run_reports/test_bc_train_dataset.md"
    model_path = "output_data/models/test_bc_policy.pt"
    train_npz = "output_data/test_bc_train_metrics.npz"
    train_report = "docs/control/run_reports/test_bc_train_report.md"

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

    assert os.path.exists(model_path)
    assert os.path.exists(train_npz)
    assert os.path.exists(train_report)
