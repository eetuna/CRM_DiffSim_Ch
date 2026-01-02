import glob
import os
import subprocess
import sys

import numpy as np


def test_replay_dataset_mode_circle(tmp_path):
    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    before = set(glob.glob(os.path.join(out_dir, "replay_mode_circle_*.npz")))
    cmd = [
        sys.executable,
        "examples/run_ilqr_circle.py",
        "--backend",
        "v1_3",
        "--replay-dataset-mode",
        "--mpc-steps",
        "5",
        "--horizon",
        "10",
        "--max-iter",
        "1",
        "--max-wall-sec",
        "120",
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
    assert os.path.getsize(latest) > 0

    data = np.load(latest, allow_pickle=True)
    tip = data["tip_xyz"]
    target = data["target_xyz"]
    assert np.isfinite(tip).all()
    assert np.isfinite(target).all()
    n = min(tip.shape[0], target.shape[0])
    errors = np.linalg.norm(tip[:n] - target[:n], axis=1)
    assert float(errors.mean()) < 50.0
