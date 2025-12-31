from pathlib import Path

import numpy as np
import pytest
import torch

from crm_diffsims.control.ilqr import ILQRConfig, run_mpc, tip_position
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_ilqr_mpc_smoke(monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    build_dir = repo_root / "build_torch_ext_dyn"
    if not any(build_dir.glob("crm_dynamics_ext*.so")):
        pytest.skip("crm_dynamics_ext not built yet; run warmup to compile extension")

    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    x0, li = build_state()
    tip0 = tip_position(x0)[0]
    target_offset = torch.tensor([0.0, 0.2, 0.0], dtype=torch.float64)

    n_steps = 1
    targets = (tip0 + target_offset).repeat(n_steps + 1, 1)

    horizon = 3
    u_warm = torch.zeros((horizon, 1, 3), dtype=torch.float64)

    cfg = ILQRConfig(
        horizon=horizon,
        max_iter=1,
        w_tip=1.0,
        w_tip_terminal=5.0,
        w_u=1e-2,
        w_du=1e-1,
        max_du=0.01,
        max_sign_flip=0.005,
    )

    def _linearize_fast(x_t, u_t, li_in, cfg_in):
        n_x = x_t.shape[-1]
        n_u = u_t.shape[-1]
        a = torch.eye(n_x, dtype=x_t.dtype, device=x_t.device)
        b = torch.zeros((n_x, n_u), dtype=x_t.dtype, device=x_t.device)
        return a, b

    monkeypatch.setattr("crm_diffsims.control.ilqr._linearize", _linearize_fast)

    x_roll, u_roll, _, stats = run_mpc(x0, targets, li, cfg_dyn, cfg, u_warm)

    assert torch.isfinite(x_roll).all().item()
    assert torch.isfinite(u_roll).all().item()

    for idx in range(u_roll.shape[0]):
        x_next, _, exit_flag, _, _, _ = ext.crm_step_forward(x_roll[idx:idx + 1], u_roll[idx:idx + 1], li, cfg_dyn)
        assert int(exit_flag.item()) == 0
        assert torch.isfinite(x_next).all().item()

    clamp_rates = [s["clamp_rate"] for s in stats]
    assert max(clamp_rates) < 1.0
