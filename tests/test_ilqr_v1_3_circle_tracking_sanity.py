import os
import sys

import numpy as np
import torch

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
python_root = os.path.join(repo_root, "python")
if python_root not in sys.path:
    sys.path.insert(0, python_root)

from crm_diffsims.control import ilqr as ilqr_mod
from crm_diffsims.control.ilqr import ILQRConfig, tip_position
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_ilqr_v1_3_circle_tracking_sanity():
    torch.manual_seed(0)
    np.random.seed(0)

    data = np.load("data/dyn_fk_ramp_circle1_hold1.npz")
    targets_np = data["tip_dyn"]
    currents_np = data["currents"]

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, _ = build_state()
    li = torch.tensor([float(data.get("insertion_length", 94.3))], dtype=torch.float64)

    n_steps = 2
    targets = torch.tensor(targets_np[: n_steps + 1], dtype=torch.float64)
    horizon = 3

    cfg = ILQRConfig(
        horizon=horizon,
        max_iter=2,
        w_tip=1.0,
        w_tip_terminal=25.0,
        w_u=5e-3,
        w_du=5e-2,
        max_du=0.008,
        max_sign_flip=0.004,
    )
    cfg.linearization_backend = "v1_3"
    cfg.linearization_dense = True

    safe_max_u = cfg.max_u
    safe_max_du = cfg.max_du
    safe_path = "docs/control/safe_bounds.json"
    if os.path.exists(safe_path):
        with open(safe_path, "r", encoding="ascii") as f:
            import json

            bounds = json.load(f)
        safe_max_u = float(bounds.get("umax_safe", cfg.max_u))
        safe_max_du = float(bounds.get("d_umax_safe", cfg.max_du))

    cfg.max_u = min(safe_max_u, 0.1)
    cfg.max_du = min(safe_max_du, 0.01)

    u_seed = torch.tensor(currents_np[:horizon], dtype=torch.float64).view(horizon, 1, 3)
    u0 = torch.zeros((1, 1, 3), dtype=torch.float64)
    ramped = []
    for t in range(u_seed.shape[0]):
        alpha = float(t + 1) / max(1, u_seed.shape[0])
        ramped.append((1.0 - alpha) * u0 + alpha * u_seed[t : t + 1])
    u_seq = torch.cat(ramped, dim=0)

    u_clamped = u_seq.clone()
    for t in range(u_clamped.shape[0]):
        u_clamped[t] = torch.clamp(u_clamped[t], -cfg.max_u, cfg.max_u)
        if t > 0:
            delta = u_clamped[t] - u_clamped[t - 1]
            delta_clamped = torch.clamp(delta, -cfg.max_du, cfg.max_du)
            sign_flip = (u_clamped[t - 1] * (u_clamped[t - 1] + delta_clamped)) < 0.0
            if sign_flip.any():
                delta_clamped = torch.where(
                    sign_flip,
                    torch.clamp(delta_clamped, -cfg.max_sign_flip, cfg.max_sign_flip),
                    delta_clamped,
                )
            u_clamped[t] = u_clamped[t - 1] + delta_clamped
    u_seq = u_clamped

    x_roll = [x0]
    errors = []
    for t_step in range(n_steps):
        alpha = min(1.0, max(0.0, t_step / max(1, 1)))
        cfg.max_u = min(safe_max_u, 0.1 + (0.3 - 0.1) * alpha)
        cfg.max_du = min(safe_max_du, 0.01 + (0.05 - 0.01) * alpha)
        horizon_targets = targets[t_step : t_step + cfg.horizon + 1]
        if horizon_targets.shape[0] < cfg.horizon + 1:
            pad = horizon_targets[-1:].repeat(cfg.horizon + 1 - horizon_targets.shape[0], 1)
            horizon_targets = torch.cat([horizon_targets, pad], dim=0)
        result = ilqr_mod.ilqr_solve(x_roll[-1], u_seq, horizon_targets, li, cfg_dyn, cfg)
        u_apply = result.u_seq[0:1]
        x_next, step_ok = ilqr_mod._step_forward(
            x_roll[-1],
            u_apply,
            li,
            cfg_dyn,
            cfg.linearization_backend,
        )
        assert step_ok
        assert torch.isfinite(x_next).all().item()
        x_roll.append(x_next)
        u_seq = torch.cat([result.u_seq[1:], result.u_seq[-1:]], dim=0)
        tip_now = tip_position(x_next)[0].detach().cpu().numpy()
        target_now = targets[t_step + 1].detach().cpu().numpy()
        errors.append(np.linalg.norm(tip_now - target_now))

    assert np.isfinite(errors).all()
    assert errors[-1] < errors[0]
