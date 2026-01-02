import os
import subprocess
import sys

import numpy as np


def test_bc_dataset_smoke():
    out_npz = "output_data/test_bc_dataset_smoke.npz"
    report_path = "docs/control/run_reports/test_bc_dataset_smoke.md"

    if os.path.exists(out_npz):
        os.remove(out_npz)
    if os.path.exists(report_path):
        os.remove(report_path)

    cmd = [
        sys.executable,
        "examples/generate_bc_dataset.py",
        "--controller",
        "replay_only_full",
        "--episodes",
        "1",
        "--episode-len",
        "6",
        "--start-idx",
        "0",
        "--end-idx",
        "12",
        "--seed",
        "0",
        "--output-npz",
        out_npz,
        "--report-path",
        report_path,
    ]
    subprocess.run(cmd, check=True)

    assert os.path.exists(out_npz)
    assert os.path.exists(report_path)

    with np.load(out_npz) as data:
        assert "x_t" in data
        assert "target_xyz" in data
        assert "u_t" in data
        assert data["x_t"].shape[0] > 0
