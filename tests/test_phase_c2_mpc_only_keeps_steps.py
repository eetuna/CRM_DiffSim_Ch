import os
import subprocess
import sys

import numpy as np


def test_phase_c2_mpc_only_keeps_steps():
    env = os.environ.copy()
    env["PYTHONPATH"] = "."
    env["TORCH_EXTENSIONS_DIR"] = "/workspaces/CRM_DiffSim_Ch/build_torch_ext_dyn"

    eval_cmd = [
        sys.executable,
        "examples/eval_policy_warmstart_mpc.py",
        "--model",
        "output_data/models/phase_c2_policy.pt",
        "--dataset",
        "data/dyn_fk_ramp_circle1_hold1.npz",
        "--start-idx",
        "0",
        "--window-len",
        "50",
        "--init-mode",
        "dataset_x0",
        "--preroll-steps",
        "5",
        "--horizon",
        "5",
        "--mpc-steps",
        "3",
        "--max-wall-sec",
        "120",
        "--seed",
        "0",
        "--output-prefix",
        "phase_c2_warmstart_mpc_only_test",
        "--report-path",
        "docs/control/run_reports/phase_c2_warmstart_mpc_only_test.md",
    ]

    proc = subprocess.run(eval_cmd, check=False, env=env, timeout=300)
    assert proc.returncode in (0, 2)

    npz_path = "output_data/phase_c2_warmstart_mpc_only_test.npz"
    assert os.path.exists(npz_path)

    with np.load(npz_path, allow_pickle=False) as data:
        keep = data["keep_mpc_only"].astype(bool)
        assert keep.sum() > 0
        err = data["error_mpc_only"]
        if keep.sum() > 0:
            assert np.isfinite(err).all()
