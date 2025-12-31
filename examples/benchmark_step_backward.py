import os
import time

import torch

from crm_diffsims.dynamics.step import build_cfg, build_state, crm_step, load_dyn_ext


def _rollout(x0, u_seq, li, cfg):
    x_t = x0
    xs = []
    etas = []
    for t in range(u_seq.shape[0]):
        x_tp1, eta, exit_flag, res_norm, u0, tau = ext.crm_step_forward(x_t, u_seq[t:t + 1], li, cfg)
        xs.append(x_t)
        etas.append(eta)
        x_t = x_tp1
    return xs, etas


def _backward_fd(xs, etas, u_seq, li, cfg):
    os.environ["CRM_DIFFSIMS_ENABLE_FD_REFERENCE"] = "1"
    for x_t, eta, u_t in zip(xs, etas, u_seq):
        ext.crm_dyn_jacobians_fd(eta, x_t, u_t.unsqueeze(0), li, cfg)
        grad = torch.ones_like(x_t)
        ext.crm_dyn_vjp_fd(grad, eta, x_t, u_t.unsqueeze(0), li, cfg)


def _backward_new(xs, etas, u_seq, li, cfg):
    for x_t, eta, u_t in zip(xs, etas, u_seq):
        ext.crm_dyn_jacobians(eta, x_t, u_t.unsqueeze(0), li, cfg)
        grad = torch.ones_like(x_t)
        ext.crm_dyn_vjp(grad, eta, x_t, u_t.unsqueeze(0), li, cfg)


def main():
    torch.manual_seed(0)
    os.environ.setdefault("CRM_DIFFSIMS_ENABLE_FD_REFERENCE", "1")

    global ext
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x0, li = build_state()
    u_seq = torch.tensor(
        [
            [[0.0, 0.0, 0.03]],
            [[0.01, -0.005, 0.035]],
            [[-0.01, 0.005, 0.032]],
            [[0.005, 0.0, 0.03]],
            [[-0.008, 0.004, 0.033]],
            [[0.0, -0.006, 0.031]],
            [[0.007, 0.003, 0.034]],
            [[-0.006, -0.003, 0.03]],
            [[0.01, 0.005, 0.035]],
            [[-0.01, -0.005, 0.032]],
        ],
        dtype=torch.float64,
    )

    xs, etas = _rollout(x0, u_seq, li, cfg)

    start = time.perf_counter()
    _backward_fd(xs, etas, u_seq, li, cfg)
    fd_time = time.perf_counter() - start

    start = time.perf_counter()
    _backward_new(xs, etas, u_seq, li, cfg)
    new_time = time.perf_counter() - start

    print(f"fd_reference_time={fd_time:.6f}s")
    print(f"v1_2_time={new_time:.6f}s")


if __name__ == "__main__":
    main()
