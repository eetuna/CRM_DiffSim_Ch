import os
import subprocess
import sys

import numpy as np


def test_generate_bc_dataset_cem_smoke():
    out_npz = "output_data/test_bc_dataset_cem_smoke.npz"
    report_path = "docs/control/run_reports/test_bc_dataset_cem_smoke.md"

    if os.path.exists(out_npz):
        os.remove(out_npz)
    if os.path.exists(report_path):
        os.remove(report_path)

    cmd = [
        sys.executable,
        "examples/generate_bc_dataset.py",
        "--controller",
        "cem_mpc",
        "--episodes",
        "1",
        "--episode-len",
        "20",
        "--start-idx",
        "0",
        "--end-idx",
        "20",
        "--horizon",
        "5",
        "--mpc-steps",
        "5",
        "--num-samples",
        "8",
        "--num-elites",
        "2",
        "--cem-iters",
        "1",
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

    with np.load(out_npz, allow_pickle=False) as data:
        assert "keep_mask" in data
        assert data["keep_mask"].sum() > 0
