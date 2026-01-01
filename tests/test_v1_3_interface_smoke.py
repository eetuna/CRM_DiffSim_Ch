import torch

from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


def test_v1_3_interface_smoke():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x0, li = build_state()
    u0 = torch.zeros((1, 1, 3), dtype=torch.float64)

    y1, cache1 = crm_step_v1_3_forward(x0, u0, li, cfg)
    y2, cache2 = crm_step_v1_3_forward(x0, u0, li, cfg)

    assert y1.shape == y2.shape
    assert y1.dtype == torch.float64
    assert torch.allclose(y1, y2)

    for key in ["eta", "solver_exit", "residual_norm", "u0", "tau"]:
        assert key in cache1
        assert torch.allclose(cache1[key], cache2[key])
