import torch
import pytest

from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3


def test_v1_3_backward_nonconvergence():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x0, li = build_state()
    u = torch.tensor([[[5.0, 5.0, 5.0]]], dtype=torch.float64, requires_grad=True)
    x0 = x0.clone().detach().requires_grad_(True)

    y = crm_step_v1_3(x0, u, li, cfg)
    grad = torch.ones_like(y)

    with pytest.raises(RuntimeError):
        torch.autograd.grad(y, (x0, u), grad_outputs=grad, retain_graph=True)
