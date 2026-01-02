import os

import numpy as np
import torch

from crm_diffsims.control.ilqr import ILQRConfig, run_mpc, tip_position
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def _build_staged_targets(full_targets, base_tip, n_steps, n_warm, n_blend):
    staged_targets = full_targets.clone()
    warm_end = min(n_warm, n_steps)
    for t in range(warm_end + 1):
        staged_targets[t] = base_tip
    blend_start = warm_end + 1
    blend_end = min(n_steps, warm_end + max(n_blend, 0))
    if n_blend > 0 and blend_start <= blend_end:
        steps_count = blend_end - blend_start + 1
        for idx, t in enumerate(range(blend_start, blend_end + 1)):
            if steps_count > 1:
                alpha = idx / (steps_count - 1)
            else:
                alpha = 1.0
            staged_targets[t] = (1.0 - alpha) * base_tip + alpha * full_targets[t]
    return staged_targets


def test_ilqr_mpc_v1_3_lemniscate_staged_smoke():
    torch.manual_seed(0)
    np.random.seed(0)

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, _ = build_state()
    li = torch.tensor([94.3], dtype=torch.float64)

    target_np = np.load("data/dyn_fk_lem1_y40_a10_L94_hold1.npz")["tip_dyn"]
    n_steps = 2
    full_targets = torch.tensor(target_np[: n_steps + 1], dtype=torch.float64)
    base_tip = tip_position(x0)[0]

    staged_targets = _build_staged_targets(full_targets, base_tip, n_steps, n_warm=0, n_blend=1)

    cfg = ILQRConfig(horizon=3, max_iter=1, max_du=0.008, max_sign_flip=0.004)
    cfg.linearization_backend = "v1_3"
    cfg.linearization_dense = True

    safe_path = "docs/control/safe_bounds.json"
    if os.path.exists(safe_path):
        with open(safe_path, "r", encoding="ascii") as f:
            import json

            bounds = json.load(f)
        cfg.max_u = float(bounds.get("umax_safe", cfg.max_u))
        cfg.max_du = float(bounds.get("d_umax_safe", cfg.max_du))

    u_warm = torch.zeros((cfg.horizon, 1, 3), dtype=torch.float64)

    x_roll, u_roll, u_final, stats = run_mpc(x0, staged_targets, li, cfg_dyn, cfg, u_warm)

    assert torch.isfinite(x_roll).all().item()
    assert torch.isfinite(u_roll).all().item()
    assert stats

    max_u = float(cfg.max_u)
    max_du = float(cfg.max_du)
    u_abs_max = u_roll.abs().max().item() if u_roll.numel() else 0.0
    assert u_abs_max <= max_u + 1e-6
    if u_roll.shape[0] > 1:
        du_abs_max = (u_roll[1:] - u_roll[:-1]).abs().max().item()
        assert du_abs_max <= max_du + 1e-6
