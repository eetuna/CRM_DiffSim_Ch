import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, crm_step, load_dyn_ext


def test_crm_step_gradcheck():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x_t, li = build_stable_state(ext, cfg)
    x_t = x_t.clone().detach().requires_grad_(True)
    u_t = torch.tensor([[[0.0, 0.0, 0.03]]], dtype=torch.float64, requires_grad=True)

    def func(x_in, u_in):
        return crm_step(x_in, u_in, li, cfg)

    assert torch.autograd.gradcheck(func, (x_t, u_t), eps=1e-6, atol=1e-4, rtol=1e-3)
