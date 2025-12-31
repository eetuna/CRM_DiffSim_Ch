import numpy as np
import torch

from crm_diffsims.control.ilqr import ILQRConfig, run_mpc, tip_position
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def test_ilqr_full_fidelity_smoke():
    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    x0, li = build_state()
    tip0 = tip_position(x0)[0]
    target_offset = torch.tensor([0.0, 0.1, 0.0], dtype=torch.float64)

    horizon = 1
    n_steps = 1
    targets = (tip0 + target_offset).repeat(n_steps + 1, 1)
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
        line_search_alphas=(1.0,),
    )

    x_roll, u_roll, _, stats = run_mpc(x0, targets, li, cfg_dyn, cfg, u_warm)

    assert torch.isfinite(x_roll).all().item()
    assert torch.isfinite(u_roll).all().item()
    assert max(s["clamp_rate"] for s in stats) < 1.0
