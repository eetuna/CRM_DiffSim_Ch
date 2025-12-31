import torch

from crm_diffsims.dynamics.step import build_cfg, build_state, crm_step, load_dyn_ext


def main():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)

    x_t, li = build_state()
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
        requires_grad=True,
    )

    tip_idx = [12, 13, 14, 24, 25, 26]
    loss = torch.zeros((), dtype=torch.float64)
    for t in range(u_seq.shape[0]):
        x_t = crm_step(x_t, u_seq[t:t + 1], li, cfg)
        loss = loss + (x_t[:, tip_idx] ** 2).sum()

    loss.backward()

    grad_norm = u_seq.grad.norm().item()
    print("crm_step rollout/backward OK")
    print(f"loss={loss.item():.6e}")
    print(f"grad_norm={grad_norm:.6e}")


if __name__ == "__main__":
    main()
