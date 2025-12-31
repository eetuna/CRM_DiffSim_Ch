from datetime import datetime
import argparse
import os
import sys
import time

import numpy as np
import torch

from crm_diffsims.control import ilqr as ilqr_mod
from crm_diffsims.control.ilqr import ILQRConfig, run_mpc, tip_index_sanity_check, tip_position
from crm_diffsims.control.step_profiler import run_rollout_profile
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def _load_targets(path):
    data = np.load(path)
    if "tip_desired" in data:
        target = data["tip_desired"]
    elif "tip_fk" in data:
        target = data["tip_fk"]
    else:
        target = data["tip_dyn"]
    currents = data.get("currents")
    dt = float(data.get("dt", 0.2))
    li = float(data.get("insertion_length", 94.3))
    return target, currents, dt, li


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast", action="store_true", help="Run a short validation rollout")
    parser.add_argument("--profile-step", action="store_true", help="Profile crm_step and save flight recorder")
    parser.add_argument("--max-wall-sec", type=float, default=None, help="Hard wall time limit in seconds")
    parser.add_argument("--save-u-seq", action="store_true", help="Save iLQR u_seq after MPC run")
    args = parser.parse_args()

    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    target_np, currents_np, dt, li_val = _load_targets("data/dyn_fk_lem1_y40_a10_L94_hold1.npz")

    n_steps = 100
    targets = torch.tensor(target_np[: n_steps + 1], dtype=torch.float64)
    x0, _ = build_state()
    li = torch.tensor([li_val], dtype=torch.float64)

    horizon = 35
    if currents_np is None:
        u_warm = torch.zeros((horizon, 1, 3), dtype=torch.float64)
    else:
        u_seed = torch.tensor(currents_np[:horizon], dtype=torch.float64).view(horizon, 1, 3)
        u_warm = u_seed.clone()

    cfg = ILQRConfig(
        horizon=horizon,
        max_iter=8,
        w_tip=1.0,
        w_tip_terminal=30.0,
        w_u=5e-3,
        w_du=6e-2,
        max_du=0.008,
        max_sign_flip=0.004,
    )

    if args.fast:
        n_steps = min(n_steps, 3)
        base_tip = tip_position(x0)[0]
        targets = base_tip.repeat(n_steps + 1, 1)
        cfg.horizon = min(cfg.horizon, 10)
        cfg.max_iter = min(cfg.max_iter, 5)
        u_warm = torch.zeros((cfg.horizon, 1, 3), dtype=torch.float64)
        def _linearize_fast(x_t, u_t, li_in, cfg_in):
            n_x = x_t.shape[-1]
            n_u = u_t.shape[-1]
            a = torch.eye(n_x, dtype=x_t.dtype, device=x_t.device)
            b = torch.zeros((n_x, n_u), dtype=x_t.dtype, device=x_t.device)
            return a, b

        ilqr_mod._linearize = _linearize_fast
        print("FAST mode: using approximate linearization for quick validation")

    tip_index_sanity_check(x0, u_warm[:1], li, cfg_dyn)

    start_time = time.time()
    x_roll, u_roll, u_final, stats = run_mpc(
        x0,
        targets,
        li,
        cfg_dyn,
        cfg,
        u_warm,
        max_wall_sec=args.max_wall_sec,
    )
    tip_roll = tip_position(x_roll)

    t = np.arange(tip_roll.shape[0]) * dt
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    if args.save_u_seq:
        out_dir = "output_data"
        os.makedirs(out_dir, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        first_iter_u = None
        if stats and stats[0].get("first_iter_u_seq") is not None:
            first_iter_u = stats[0]["first_iter_u_seq"].detach().cpu().numpy()
        np.savez(
            f"{out_dir}/ilqr_lemniscate_u_seq_{stamp}.npz",
            u_seq=u_final.detach().cpu().numpy(),
            first_iter_u_seq=first_iter_u,
        )

    if args.max_wall_sec is not None and time.time() - start_time > args.max_wall_sec:
        out_dir = "output_data"
        os.makedirs(out_dir, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        np.savez(
            f"{out_dir}/ilqr_lemniscate_partial_{stamp}.npz",
            x_roll=x_roll.detach().cpu().numpy(),
            u_seq=u_final.detach().cpu().numpy(),
        )
        if args.profile_step and u_roll.shape[0] > 0:
            run_rollout_profile(x0, u_roll, li, cfg_dyn, label="lemniscate_profile_partial", save_dir=out_dir)
        print("Max wall time exceeded, saved partial rollout.")
        sys.exit(2)

    if any(s.get("nan_detected") for s in stats):
        out_dir = "output_data"
        os.makedirs(out_dir, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        np.savez(
            f"{out_dir}/ilqr_lemniscate_nan_{stamp}.npz",
            x_roll=x_roll.detach().cpu().numpy(),
            u_seq=u_final.detach().cpu().numpy(),
            iter_index=0,
            fail_step=-1,
        )
        print("NaN detected in iLQR; saved diagnostics.")
        sys.exit(1)

    if args.fast:
        print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
        return

    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError as exc:
        raise SystemExit("matplotlib is required for plotting") from exc

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(targets[:, 0].cpu().numpy(), targets[:, 1].cpu().numpy(), targets[:, 2].cpu().numpy(), label="target")
    ax.plot(tip_roll[:, 0].cpu().numpy(), tip_roll[:, 1].cpu().numpy(), tip_roll[:, 2].cpu().numpy(), label="ilqr")
    ax.set_title("Tip trajectory (lemniscate)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend()

    fig2, ax2 = plt.subplots(figsize=(8, 4))
    u_np = u_roll[:, 0].cpu().numpy()
    ax2.plot(t[: u_np.shape[0]], u_np[:, 0], label="i1")
    ax2.plot(t[: u_np.shape[0]], u_np[:, 1], label="i2")
    ax2.plot(t[: u_np.shape[0]], u_np[:, 2], label="i3")
    ax2.set_title("Currents (lemniscate)")
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("current")
    ax2.legend()

    errors = np.linalg.norm(tip_roll[: targets.shape[0]].cpu().numpy() - targets.cpu().numpy(), axis=1)
    fig3, ax3 = plt.subplots(figsize=(8, 4))
    ax3.plot(t[: errors.shape[0]], errors, label="tip error")
    ax3.set_title("Tip tracking error (lemniscate)")
    ax3.set_xlabel("time [s]")
    ax3.set_ylabel("error")
    ax3.legend()

    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    fig.savefig(f"{out_dir}/ilqr_lemniscate_tip_{stamp}.png", dpi=150)
    fig2.savefig(f"{out_dir}/ilqr_lemniscate_currents_{stamp}.png", dpi=150)
    fig3.savefig(f"{out_dir}/ilqr_lemniscate_error_{stamp}.png", dpi=150)
    print(f"Saved plots to {out_dir}/ilqr_lemniscate_*_{stamp}.png")
    print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
    for step in range(1, errors.shape[0]):
        mean_err = errors[: step + 1].mean()
        final_err = errors[step]
        print(f"MPC step {step}: mean_tip_error={mean_err:.3e} final_tip_error={final_err:.3e}")
    plt.show()

    if args.profile_step:
        run_rollout_profile(x0, u_roll, li, cfg_dyn, label="lemniscate_profile", save_dir="output_data")


if __name__ == "__main__":
    main()
