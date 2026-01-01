import json
import os
from pathlib import Path

import numpy as np
import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext


def _load_safe_bounds():
    path = Path("docs/control/safe_bounds.json")
    if not path.exists():
        return 0.1, 0.01
    with path.open("r", encoding="ascii") as f:
        data = json.load(f)
    return float(data.get("umax_safe", 0.1)), float(data.get("d_umax_safe", 0.01))


def _fd_col_residual(ext, eta, x, u, li, cfg, col, kind, eps):
    if kind == "eta":
        plus = eta.clone()
        minus = eta.clone()
        plus[0, col] += eps
        minus[0, col] -= eps
        r_plus = ext.crm_dyn_residual(plus, x, u, li, cfg)[0]
        r_minus = ext.crm_dyn_residual(minus, x, u, li, cfg)[0]
    elif kind == "x":
        plus = x.clone()
        minus = x.clone()
        plus[0, col] += eps
        minus[0, col] -= eps
        r_plus = ext.crm_dyn_residual(eta, plus, u, li, cfg)[0]
        r_minus = ext.crm_dyn_residual(eta, minus, u, li, cfg)[0]
    elif kind == "u":
        plus = u.clone()
        minus = u.clone()
        plus[0].reshape(-1)[col] += eps
        minus[0].reshape(-1)[col] -= eps
        r_plus = ext.crm_dyn_residual(eta, x, plus, li, cfg)[0]
        r_minus = ext.crm_dyn_residual(eta, x, minus, li, cfg)[0]
    else:
        raise ValueError(f"unknown kind: {kind}")
    return (r_plus - r_minus) / (2.0 * eps)


def _stable_fd_col(ext, eta, x, u, li, cfg, col, kind, base, eps_list):
    cols = []
    for eps in eps_list:
        step = eps * (1.0 + 0.01 * abs(base))
        cols.append(_fd_col_residual(ext, eta, x, u, li, cfg, col, kind, step))
    cols = torch.stack(cols)
    mean = cols.mean(dim=0)
    rel_spread = (cols - mean).abs().max() / (mean.abs().max() + 1e-12)
    assert rel_spread.item() < 0.2
    return mean


def test_v1_3_blocks_vs_fd_local():
    torch.manual_seed(0)
    np.random.seed(0)
    os.environ["CRM_DIFFSIMS_ENABLE_FD_REFERENCE"] = "1"

    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x0, li = build_stable_state(ext, cfg)
    umax_safe, d_umax_safe = _load_safe_bounds()

    u = torch.zeros((1, 1, 3), dtype=torch.float64)
    u[0, 0, 2] = 0.5 * d_umax_safe

    y, eta, solver_exit, residual_norm, u0, tau = ext.crm_step_forward(x0, u, li, cfg)
    assert int(solver_exit.item()) == 0
    assert float(residual_norm.item()) <= cfg.residual_threshold

    A_step, B_step, C_step = ext.crm_dyn_jacobians(eta, x0, u, li, cfg)
    J_eta = A_step[0]
    J_u = B_step[0]
    J_x = C_step[0]

    eps_list = [1e-6, 5e-6, 1e-5]

    # Sample a few columns to keep runtime bounded.
    eta_cols = [0, 2, 5]
    for col in eta_cols:
        base = float(eta[0, col].item())
        fd = _stable_fd_col(ext, eta, x0, u, li, cfg, col, "eta", base, eps_list)
        diff = (J_eta[:, col] - fd).abs().max().item()
        ref = fd.abs().max().item() + 1e-12
        assert diff / ref < 1e-2

    x_cols = [0, 3, 10]
    for col in x_cols:
        base = float(x0[0, col].item())
        fd = _stable_fd_col(ext, eta, x0, u, li, cfg, col, "x", base, eps_list)
        diff = (J_x[:, col] - fd).abs().max().item()
        ref = fd.abs().max().item() + 1e-12
        assert diff / ref < 1e-2

    u_cols = [0, 2]
    for col in u_cols:
        base = float(u[0].reshape(-1)[col].item())
        fd = _stable_fd_col(ext, eta, x0, u, li, cfg, col, "u", base, eps_list)
        diff = (J_u[:, col] - fd).abs().max().item()
        ref = fd.abs().max().item() + 1e-12
        assert diff / ref < 1e-2

    # Local FD on h(eta, x, u) via vjp FD reference.
    v = torch.randn_like(y)
    vjp_eta, vjp_u, vjp_x = ext.crm_dyn_vjp(v, eta, x0, u, li, cfg)
    vjp_eta_fd, vjp_u_fd, vjp_x_fd = ext.crm_dyn_vjp_fd(v, eta, x0, u, li, cfg)

    assert (vjp_eta - vjp_eta_fd).abs().max().item() < 1e-6
    assert (vjp_u - vjp_u_fd).abs().max().item() < 1e-6
    assert (vjp_x - vjp_x_fd).abs().max().item() < 1e-6
