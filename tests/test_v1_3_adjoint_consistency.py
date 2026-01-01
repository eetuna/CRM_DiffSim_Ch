import json
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


def test_v1_3_adjoint_consistency():
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

    vjp_eta_direct, vjp_u_direct, vjp_x_direct = ext.crm_dyn_vjp(v, eta, x0, u, li, cfg)
    A_step, B_step, C_step = ext.crm_dyn_jacobians(eta, x0, u, li, cfg)

    J_eta = A_step[0]
    J_u = B_step[0]
    J_x = C_step[0]

    rhs = vjp_eta_direct[0].unsqueeze(1)
    lambda_vec = torch.linalg.solve(J_eta.transpose(0, 1), rhs).squeeze(1)

    residual = J_eta.transpose(0, 1).matmul(lambda_vec) - rhs.squeeze(1)
    denom = rhs.squeeze(1).norm() + 1e-12
    rel_res = residual.norm() / denom
    assert rel_res.item() < 1e-6

    vjp_x_ref = vjp_x_direct[0] - J_x.transpose(0, 1).matmul(lambda_vec)
    vjp_u_ref = vjp_u_direct[0].reshape(-1) - J_u.transpose(0, 1).matmul(lambda_vec)

    vjp_x, vjp_u = ext.crm_dyn_vjp_v1_3(v, eta, x0, u, li, cfg)
    vjp_u_flat = vjp_u[0].reshape(-1)

    abs_x = (vjp_x[0] - vjp_x_ref).abs().max().item()
    abs_u = (vjp_u_flat - vjp_u_ref).abs().max().item()
    rel_x = abs_x / (vjp_x_ref.abs().max().item() + 1e-12)
    rel_u = abs_u / (vjp_u_ref.abs().max().item() + 1e-12)

    assert abs_x < 1e-7 or rel_x < 1e-6
    assert abs_u < 1e-7 or rel_u < 1e-6
