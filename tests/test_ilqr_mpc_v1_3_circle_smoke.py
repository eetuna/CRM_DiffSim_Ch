import time

import numpy as np
import torch

from crm_diffsims.control.ilqr import ILQRConfig, run_mpc
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_ilqr_mpc_v1_3_circle_smoke():
    torch.manual_seed(0)
    np.random.seed(0)

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, _ = build_state()
    li = torch.tensor([94.3], dtype=torch.float64)

    target_np = np.load("data/dyn_fk_ramp_circle1_hold1.npz")["tip_dyn"]
    n_steps = 3
    base_targets = torch.tensor(target_np[: n_steps + 1], dtype=torch.float64)

    cfg = ILQRConfig(horizon=5, max_iter=1, max_du=0.008, max_sign_flip=0.004)
    cfg.linearization_backend = "v1_3"
    cfg.linearization_dense = True

    u_warm = torch.zeros((cfg.horizon, 1, 3), dtype=torch.float64)

    horizon_targets = base_targets
    if horizon_targets.shape[0] < cfg.horizon + 1:
        pad = horizon_targets[-1:].repeat(cfg.horizon + 1 - horizon_targets.shape[0], 1)
        horizon_targets = torch.cat([horizon_targets, pad], dim=0)

    start = time.time()
    x_roll, u_roll, u_final, stats = run_mpc(x0, base_targets, li, cfg_dyn, cfg, u_warm)
    runtime = time.time() - start
    print(f"circle v1_3 smoke runtime_sec={runtime:.2f}")

    assert torch.isfinite(x_roll).all().item()
    assert all(not s["nan_detected"] for s in stats)
    assert any("clamp_rate" in s for s in stats)

    costs = [s["cost"] for s in stats if "cost" in s]
    assert all(np.isfinite(costs))
    if len(costs) > 1:
        assert costs[-1] <= costs[0] * 1.1
