import argparse
import time

import numpy as np
import torch

from crm_diffsims.control.cem_mpc import CEMConfig, cem_mpc, tip_position as cem_tip
from crm_diffsims.control.ilqr import ILQRConfig, run_mpc, tip_position as ilqr_tip
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def _load_targets(path):
    data = np.load(path)
    if "tip_desired" in data:
        target = data["tip_desired"]
    elif "tip_fk" in data:
        target = data["tip_fk"]
    else:
        target = data["tip_dyn"]
    dt = float(data.get("dt", 0.2))
    return target, dt


def _metrics(tip_roll, targets):
    errors = np.linalg.norm(tip_roll[: targets.shape[0]] - targets, axis=1)
    return float(errors.mean()), float(errors[-1])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="data/dyn_fk_ramp_circle1_hold1.npz")
    parser.add_argument("--steps", type=int, default=20)
    args = parser.parse_args()

    torch.manual_seed(0)
    np.random.seed(0)

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, _ = build_state()
    li = torch.tensor([94.3], dtype=torch.float64)

    target_np, dt = _load_targets(args.dataset)
    targets = torch.tensor(target_np[: args.steps + 1], dtype=torch.float64)

    cem_cfg = CEMConfig(horizon=10, mpc_steps=args.steps, num_samples=32, num_elites=8, cem_iters=2)
    start = time.time()
    cem_result = cem_mpc(x0, targets, li, cfg_dyn, cem_cfg)
    cem_time = time.time() - start
    cem_tip_roll = cem_tip(cem_result.x_roll).detach().cpu().numpy()
    cem_mean, cem_final = _metrics(cem_tip_roll, targets.cpu().numpy())

    ilqr_cfg = ILQRConfig(horizon=10, max_iter=5, max_du=0.008, max_sign_flip=0.004)
    ilqr_cfg.linearization_backend = "v1_3"
    ilqr_cfg.linearization_dense = True
    u_warm = torch.zeros((ilqr_cfg.horizon, 1, 3), dtype=torch.float64)
    start = time.time()
    x_roll, u_roll, u_final, stats = run_mpc(x0, targets, li, cfg_dyn, ilqr_cfg, u_warm)
    ilqr_time = time.time() - start
    ilqr_tip_roll = ilqr_tip(x_roll).detach().cpu().numpy()
    ilqr_mean, ilqr_final = _metrics(ilqr_tip_roll, targets.cpu().numpy())

    print("CEM vs iLQR v1.3")
    print(f"CEM:   mean_err={cem_mean:.3e} final_err={cem_final:.3e} runtime_sec={cem_time:.2f}")
    print(f"iLQR:  mean_err={ilqr_mean:.3e} final_err={ilqr_final:.3e} runtime_sec={ilqr_time:.2f}")


if __name__ == "__main__":
    main()
