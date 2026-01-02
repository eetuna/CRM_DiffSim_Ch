import os
import subprocess
import sys


def test_scoreboard_smoke():
    env = os.environ.copy()
    env["PYTHONPATH"] = "."
    env["TORCH_EXTENSIONS_DIR"] = "/workspaces/CRM_DiffSim_Ch/build_torch_ext_dyn"

    model_path = "output_data/models/phase_c2_policy.pt"
    modes = ["replay_only"]
    cmd = [
        sys.executable,
        "examples/scoreboard_eval.py",
        "--dataset",
        "data/dyn_fk_ramp_circle1_hold1.npz",
        "--start-idx",
        "0",
        "--window-len",
        "20",
        "--max-wall-sec",
        "60",
    ]
    if os.path.exists(model_path):
        modes.append("policy_only")
        cmd.extend(["--model", model_path])
    cmd.extend(["--modes", ",".join(modes)])

    result = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
    )

    npz_paths = []
    report_paths = []
    for line in result.stdout.splitlines():
        if line.startswith("artifact_npz="):
            npz_paths.append(line.split("=", 1)[1].strip())
        if line.startswith("report_path="):
            report_paths.append(line.split("=", 1)[1].strip())

    assert npz_paths, "no NPZ artifacts reported"
    assert report_paths, "no report artifacts reported"
    for path in npz_paths:
        assert os.path.exists(path)
    for path in report_paths:
        assert os.path.exists(path)
