import time

import numpy as np
import torch

from crm_diffsims.control import ilqr as ilqr_mod
from crm_diffsims.dynamics.step import build_cfg, build_state, crm_step, load_dyn_ext


def test_linearization_cost():
    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    x0, li = build_state()
    u0 = torch.zeros((1, 1, 3), dtype=torch.float64)

    t0 = time.perf_counter()
    _ = crm_step(x0, u0, li, cfg_dyn)
    step_time = time.perf_counter() - t0

    horizon = 1
    x_t = x0
    lin_time = 0.0
    for _ in range(horizon):
        t1 = time.perf_counter()
        _ = ilqr_mod._linearize(x_t, u0, li, cfg_dyn)
        lin_time += time.perf_counter() - t1
        x_t = crm_step(x_t, u0, li, cfg_dyn)

    ratio = lin_time / max(step_time * horizon, 1e-9)
    print(f"step_time={step_time:.4f}s lin_time={lin_time:.4f}s ratio={ratio:.2f}")
    assert np.isfinite(ratio)
