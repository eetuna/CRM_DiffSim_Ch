import os

import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext


def test_crm_step_C_step_matches_fd():
    os.environ["CRM_DIFFSIMS_ENABLE_FD_REFERENCE"] = "1"
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x_t, li = build_stable_state(ext, cfg)

    max_abs = 0.0
    max_rel = 0.0

    for _ in range(3):
        u_base = torch.tensor([[[0.0, 0.0, 0.03]]], dtype=torch.float64)
        u_t = u_base + 0.002 * torch.randn_like(u_base)

        x_tp1, eta, exit_flag, res_norm, u0, tau = ext.crm_step_forward(x_t, u_t, li, cfg)
        assert int(exit_flag.item()) == 0

        C_new = ext.crm_dyn_jacobians(eta, x_t, u_t, li, cfg)[2]
        C_fd = ext.crm_dyn_jacobians_fd(eta, x_t, u_t, li, cfg)[2]

        diff = torch.abs(C_new - C_fd)
        rel = diff / torch.clamp(torch.abs(C_fd), min=1e-8)
        max_abs = max(max_abs, diff.max().item())
        max_rel = max(max_rel, rel.max().item())

    print(f"max_abs={max_abs:.3e}")
    print(f"max_rel={max_rel:.3e}")

    # C_step_new uses scale-aware steps; allow small deviation vs fixed-step FD.
    assert max_abs <= 1e-4
    assert max_rel <= 2e-2
