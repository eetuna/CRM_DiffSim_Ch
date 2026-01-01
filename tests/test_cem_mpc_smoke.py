import time

import numpy as np
import torch

from crm_diffsims.control.cem_mpc import CEMConfig, cem_mpc
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_cem_mpc_smoke():
    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    x0, li = build_state()
    targets = torch.zeros((6, 3), dtype=torch.float64)
    targets[:] = x0[0, 24:27]

    cfg = CEMConfig(
        horizon=5,
        mpc_steps=2,
        num_samples=16,
        num_elites=4,
        cem_iters=1,
        u_std_init=0.005,
        seed=0,
    )

    start = time.time()
    result = cem_mpc(x0, targets, li, cfg_dyn, cfg)
    elapsed = time.time() - start

    assert torch.isfinite(result.x_roll).all().item()
    assert torch.isfinite(result.u_roll).all().item()

    umax = cfg.umax_safe + 1e-6
    d_umax = cfg.d_umax_safe + 1e-6
    u = result.u_roll[:, 0]
    assert torch.max(torch.abs(u)).item() <= umax
    if u.shape[0] > 1:
        du = u[1:] - u[:-1]
        assert torch.max(torch.abs(du)).item() <= d_umax

    assert elapsed < 60.0
