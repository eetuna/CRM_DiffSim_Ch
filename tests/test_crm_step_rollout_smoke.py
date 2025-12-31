import time

import torch

from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_crm_step_rollout_smoke():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x_t, li = build_state()
    u_seq = torch.tensor(
        [
            [[0.0, 0.0, 0.03]],
            [[0.01, -0.005, 0.035]],
            [[-0.01, 0.005, 0.032]],
            [[0.005, 0.0, 0.03]],
            [[-0.008, 0.004, 0.033]],
            [[0.0, -0.006, 0.031]],
            [[0.007, 0.003, 0.034]],
            [[-0.006, -0.003, 0.03]],
            [[0.01, 0.005, 0.035]],
            [[-0.01, -0.005, 0.032]],
        ],
        dtype=torch.float64,
    )

    max_res = 0.0
    start = time.perf_counter()
    for t in range(u_seq.shape[0]):
        x_tp1, eta, exit_flag, res_norm, u0, tau = ext.crm_step_forward(x_t, u_seq[t:t + 1], li, cfg)
        assert torch.isfinite(x_tp1).all()
        assert int(exit_flag.item()) == 0
        max_res = max(max_res, float(res_norm.item()))
        x_t = x_tp1
    end = time.perf_counter()

    print(f"runtime_per_step={(end - start) / u_seq.shape[0]:.6f}s")
    print(f"max_residual_norm={max_res:.6e}")
    assert max_res <= cfg.residual_threshold
