import argparse
import os
from datetime import datetime

import numpy as np
import torch

from crm_diffsims.control.cem_mpc import CEMConfig, cem_mpc, tip_position
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
    li = float(data.get("insertion_length", 94.3))
    if np.isnan(target).any():
        finite_rows = np.where(np.isfinite(target).all(axis=1))[0]
        if finite_rows.size > 0:
            first_idx = int(finite_rows[0])
            target[:first_idx] = target[first_idx]
        for i in range(1, target.shape[0]):
            nan_mask = np.isnan(target[i])
            if nan_mask.any():
                target[i, nan_mask] = target[i - 1, nan_mask]
    return target, dt, li


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--plot", action="store_true", help="Save plots to output_data/")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    target_np, dt, li_val = _load_targets("data/dyn_fk_lem1_y40_a10_L94_hold1.npz")
    targets = torch.tensor(target_np[:101], dtype=torch.float64)
    x0, _ = build_state()
    li = torch.tensor([li_val], dtype=torch.float64)

    cfg = CEMConfig(
        horizon=10,
        mpc_steps=10,
        num_samples=16,
        num_elites=4,
        cem_iters=1,
        u_std_init=0.001,
        seed=args.seed,
    )

    result = cem_mpc(x0, targets, li, cfg_dyn, cfg)
    tip_roll = tip_position(result.x_roll)

    t = np.arange(tip_roll.shape[0]) * dt
    min_len = min(tip_roll.shape[0], targets.shape[0])
    errors = np.linalg.norm(
        tip_roll[:min_len].cpu().numpy() - targets[:min_len].cpu().numpy(),
        axis=1,
    )

    for idx, stat in enumerate(result.stats, start=1):
        print(
            f"MPC step {idx}: mean_tip_error={stat['mean_tip_error']:.3e} "
            f"final_tip_error={stat['final_tip_error']:.3e} clamp_rate_u={stat['clamp_rate_u']:.2%} "
            f"clamp_rate_du={stat['clamp_rate_du']:.2%} u0={stat['u0']:.3f}"
        )
        if stat.get("nan_detected"):
            print("NaN detected in rollout; stopping early.")
            break

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    if args.plot:
        try:
            from crm_diffsims.control.plot_utils import (
                save_currents_plot,
                save_error_plot,
                save_tip_trajectory_plot,
            )
            save_tip_trajectory_plot(
                tip_roll[:min_len].cpu().numpy(),
                targets[:min_len].cpu().numpy(),
                title="Tip trajectory (CEM lemniscate)",
                out_path=f"{out_dir}/cem_lemniscate_tip_{stamp}",
            )
            save_currents_plot(
                result.u_roll[:, 0].cpu().numpy(),
                title="Currents (CEM lemniscate)",
                out_path=f"{out_dir}/cem_lemniscate_currents_{stamp}",
            )
            save_error_plot(
                errors,
                title="Tip tracking error (CEM lemniscate)",
                out_path=f"{out_dir}/cem_lemniscate_error_{stamp}",
            )
            print(f"Saved plots to {out_dir}/cem_lemniscate_*_{stamp}.png/.pdf")
        except RuntimeError:
            print("matplotlib not installed; skipping plots")


if __name__ == "__main__":
    main()
