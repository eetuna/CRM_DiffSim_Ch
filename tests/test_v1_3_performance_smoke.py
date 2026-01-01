import time

import numpy as np
import torch

from crm_diffsims.dynamics.linearize_v1_3 import linearize_v1_3
from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext
from crm_diffsims.dynamics.step import crm_step


def _time_once(fn, repeats):
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return float(np.mean(times))


def test_v1_3_performance_smoke():
    torch.manual_seed(0)
    np.random.seed(0)

    ext = load_dyn_ext()
    cfg = build_cfg(ext)
    x0, li = build_stable_state(ext, cfg)

    u = torch.zeros((1, 1, 3), dtype=torch.float64)
    u[0, 0, 2] = 0.005

    def autograd_jac():
        x_t = x0.detach().requires_grad_(True)
        u_t = u.detach().requires_grad_(True)

        def _dyn_x(x_in):
            return crm_step(x_in, u_t, li, cfg)

        def _dyn_u(u_in):
            return crm_step(x_t, u_in, li, cfg)

        try:
            a = torch.autograd.functional.jacobian(_dyn_x, x_t, create_graph=False, vectorize=True)
            b = torch.autograd.functional.jacobian(_dyn_u, u_t, create_graph=False, vectorize=True)
        except RuntimeError:
            a = torch.autograd.functional.jacobian(_dyn_x, x_t, create_graph=False, vectorize=False)
            b = torch.autograd.functional.jacobian(_dyn_u, u_t, create_graph=False, vectorize=False)
        return a, b

    def v1_3_op_apply():
        a_op, b_op = linearize_v1_3(x0, u, li, cfg, return_dense=False)
        dx = torch.randn_like(x0) * 1e-3
        du = torch.randn_like(u) * 1e-3
        _ = a_op.apply(dx)
        _ = b_op.apply(du)

    autograd_jac()
    v1_3_op_apply()

    auto_mean = _time_once(autograd_jac, repeats=1)
    v1_mean = _time_once(v1_3_op_apply, repeats=1)

    assert v1_mean < auto_mean / 50.0
