import os

import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext


def test_crm_step_A_step_matches_fd():
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

        A_analytic = ext.crm_dyn_jacobians(eta, x_t, u_t, li, cfg)[0]
        A_fd = ext.crm_dyn_jacobians_fd(eta, x_t, u_t, li, cfg)[0]

        diff = torch.abs(A_analytic - A_fd)
        rel = diff / torch.clamp(torch.abs(A_fd), min=1e-8)
        max_abs = max(max_abs, diff.max().item())
        max_rel = max(max_rel, rel.max().item())

    print(f"max_abs={max_abs:.3e}")
    print(f"max_rel={max_rel:.3e}")

    # A_step analytic uses the same residual evaluations with scale-aware steps; expect tight match.
    assert max_abs <= 1e-6
    assert max_rel <= 1e-5
