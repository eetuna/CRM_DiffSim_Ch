import os
import subprocess
import sys


def test_forward_stability_sweep_smoke():
    cmd = [
        sys.executable,
        "examples/sweep_forward_stability.py",
        "--umax_list",
        "0.05,0.1",
        "--d_umax_list",
        "0.01",
        "--T",
        "10",
        "--seed",
        "0",
    ]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert os.path.exists("output_data/forward_stability_sweep.csv")
