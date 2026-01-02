import glob
import os
import subprocess
import sys

import numpy as np


def test_full_dataset_replay_length():
    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    before = set(glob.glob(os.path.join(out_dir, "replay_mode_circle_*.npz")))
    cmd = [
        sys.executable,
        "examples/run_ilqr_circle.py",
        "--backend",
        "v1_3",
        "--mode",
        "replay_only_full",
        "--full-dataset",
        "--start-idx",
        "0",
        "--end-idx",
        "50",
        "--max-wall-sec",
        "300",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr)
    assert result.returncode == 0

    after = set(glob.glob(os.path.join(out_dir, "replay_mode_circle_*.npz")))
    new_files = sorted(after - before)
    assert new_files
    latest = new_files[-1]
    data = np.load(latest, allow_pickle=False)
    assert int(data["N_total"]) == 50
    assert data["dataset_currents_selected"].shape[0] == 50
    assert data["dataset_targets_selected"].shape[0] == 50
