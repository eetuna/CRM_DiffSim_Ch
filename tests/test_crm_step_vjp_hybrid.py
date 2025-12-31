import os

import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext


def test_crm_step_vjp_hybrid_matches_fd():
    os.environ["CRM_DIFFSIMS_ENABLE_FD_REFERENCE"] = "1"
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x_t, li = build_stable_state(ext, cfg)

    for _ in range(3):
        u_base = torch.tensor([[[0.0, 0.0, 0.03]]], dtype=torch.float64)
        u_t = u_base + 0.002 * torch.randn_like(u_base)

        x_tp1, eta, exit_flag, res_norm, u0, tau = ext.crm_step_forward(x_t, u_t, li, cfg)
        assert int(exit_flag.item()) == 0

        grad = torch.randn_like(x_tp1)
        vjp_hybrid = ext.crm_dyn_vjp(grad, eta, x_t, u_t, li, cfg)
        vjp_fd = ext.crm_dyn_vjp_fd(grad, eta, x_t, u_t, li, cfg)

        vjp_eta_h = vjp_hybrid[0]
        vjp_eta_fd = vjp_fd[0]

        max_abs = torch.max(torch.abs(vjp_eta_h - vjp_eta_fd)).item()
        assert max_abs <= 1e-6
