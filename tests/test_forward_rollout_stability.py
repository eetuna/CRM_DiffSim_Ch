import numpy as np
import torch

from crm_diffsims.control.step_profiler import run_rollout_profile
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_forward_rollout_stability(capsys):
    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    x0, li = build_state()
    horizon = 20
    ramp = torch.linspace(0.0, 0.01, steps=horizon, dtype=torch.float64)
    u_seq = torch.zeros((horizon, 1, 3), dtype=torch.float64)
    u_seq[:, 0, 2] = ramp

    result = run_rollout_profile(x0, u_seq, li, cfg_dyn, label="forward_rollout", save_dir="output_data")

    captured = capsys.readouterr()
    assert "Unbounded" not in captured.out
    assert result.fail_index is None
    assert result.nan_steps == []
    assert result.mean_step_time < 5.0
    print(f"mean_step_time={result.mean_step_time:.4f}s")
