import csv
import os
import subprocess
import sys


def test_forward_stability_capture_smoke():
    cmd = [
        sys.executable,
        "examples/sweep_forward_stability_subprocess.py",
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

    csv_path = "output_data/forward_stability_sweep_captured.csv"
    assert os.path.exists(csv_path)

    with open(csv_path, "r", encoding="ascii") as f:
        rows = list(csv.DictReader(f))

    assert len(rows) == 2
