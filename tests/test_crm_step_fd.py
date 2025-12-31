import torch

from crm_diffsims.dynamics.step import build_cfg, build_stable_state, crm_step, load_dyn_ext


def _loss_from_state(x_tp1):
    idx = [12, 13, 14, 24, 25, 26]
    return x_tp1[:, idx].sum()


def test_crm_step_fd():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x_t, li = build_stable_state(ext, cfg)
    x_t = x_t.clone().detach()

    points = [
        [0.0, 0.0, 0.03],
        [0.01, -0.005, 0.035],
        [-0.01, 0.005, 0.032],
    ]

    for p in points:
        u_t = torch.tensor([[p]], dtype=torch.float64, requires_grad=True)
        x_in = x_t.clone().detach().requires_grad_(True)

        x_tp1 = crm_step(x_in, u_t, li, cfg)
        loss = _loss_from_state(x_tp1)
        loss.backward()
        grad_u = u_t.grad.detach().clone().view(-1)

        fd = torch.zeros_like(grad_u)
        for i in range(3):
            h = 1e-6
            for attempt in range(2):
                u_plus = u_t.detach().clone()
                u_minus = u_t.detach().clone()
                u_plus[0, 0, i] += h
                u_minus[0, 0, i] -= h

                x_plus = crm_step(x_t, u_plus, li, cfg)
                x_minus = crm_step(x_t, u_minus, li, cfg)
                if torch.isfinite(x_plus).all() and torch.isfinite(x_minus).all():
                    fd[i] = (_loss_from_state(x_plus) - _loss_from_state(x_minus)) / (2.0 * h)
                    break
                h = 1e-5

        abs_err = torch.max(torch.abs(grad_u - fd)).item()
        rel_err = torch.max(torch.abs(grad_u - fd) / torch.clamp(torch.abs(fd), min=1e-8)).item()

        assert abs_err <= 1e-4
        assert rel_err <= 1e-3
