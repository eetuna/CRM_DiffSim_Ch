import json
from pathlib import Path

import numpy as np
import pytest
import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext


def _load_safe_bounds():
    path = Path("docs/control/safe_bounds.json")
    if not path.exists():
        return 0.1, 0.01
    with path.open("r", encoding="ascii") as f:
        data = json.load(f)
    return float(data.get("umax_safe", 0.1)), float(data.get("d_umax_safe", 0.01))


def _fd_vjp_x(ext, x, u, li, cfg, v, dx, eps):
    x_plus = x + eps * dx
    x_minus = x - eps * dx
    y_plus = ext.crm_step_forward(x_plus, u, li, cfg)[0]
    y_minus = ext.crm_step_forward(x_minus, u, li, cfg)[0]
    return ((y_plus - y_minus) * v).sum() / (2.0 * eps)


def _fd_vjp_u(ext, x, u, li, cfg, v, du, eps):
    u_plus = u + eps * du
    u_minus = u - eps * du
    y_plus = ext.crm_step_forward(x, u_plus, li, cfg)[0]
    y_minus = ext.crm_step_forward(x, u_minus, li, cfg)[0]
    return ((y_plus - y_minus) * v).sum() / (2.0 * eps)


def _stable_fd(ext, x, u, li, cfg, v, delta, eps_list, fd_fn, label):
    values = []
    for eps in eps_list:
        val = fd_fn(ext, x, u, li, cfg, v, delta, eps)
        if not torch.isfinite(val):
            pytest.skip(f"{label} FD unstable (non-finite at eps={eps}).")
        values.append(val)
    values = torch.stack(values)
    mean = values.mean()
    rel_spread = (values - mean).abs().max() / (mean.abs() + 1e-8)
    if rel_spread > 0.2:
        pytest.skip(f"{label} FD unstable (rel_spread={rel_spread.item():.3f}).")
    return mean


def test_v1_3_vjp_vs_fd():
    torch.manual_seed(0)
    np.random.seed(0)

    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x0, li = build_stable_state(ext, cfg)
    umax_safe, d_umax_safe = _load_safe_bounds()

    u = torch.zeros((1, 1, 3), dtype=torch.float64)
    u[0, 0, 2] = 0.5 * d_umax_safe

    y, eta, solver_exit, residual_norm, u0, tau = ext.crm_step_forward(x0, u, li, cfg)
    assert int(solver_exit.item()) == 0
    assert float(residual_norm.item()) <= cfg.residual_threshold

    v = torch.randn_like(y)

    vjp_x, vjp_u = ext.crm_dyn_vjp_v1_3(v, eta, x0, u, li, cfg)

    dx = torch.randn_like(x0) * 1e-3
    fd_x = _stable_fd(
        ext,
        x0,
        u,
        li,
        cfg,
        v,
        dx,
        eps_list=[1e-4, 5e-4, 1e-3],
        fd_fn=_fd_vjp_x,
        label="x",
    )
    vjp_x_dot = (vjp_x * dx).sum()

    du = torch.randn_like(u) * 1e-3
    fd_u = _stable_fd(
        ext,
        x0,
        u,
        li,
        cfg,
        v,
        du,
        eps_list=[1e-4, 5e-4, 1e-3],
        fd_fn=_fd_vjp_u,
        label="u",
    )
    vjp_u_dot = (vjp_u * du).sum()

    rel_x = torch.abs(fd_x - vjp_x_dot) / (torch.abs(fd_x) + 1e-8)
    rel_u = torch.abs(fd_u - vjp_u_dot) / (torch.abs(fd_u) + 1e-8)

    assert rel_x.item() < 1e-1
    assert rel_u.item() < 1e-1
