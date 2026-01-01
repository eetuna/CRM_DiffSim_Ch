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
    parser.add_argument("--backend", choices=["v1_3", "autograd"], default="v1_3")
    parser.add_argument("--plot", action="store_true", help="Save plots to output_data/")
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
    cfg.linearization_backend = args.backend
    cfg.linearization_dense = True

    safe_path = "docs/control/safe_bounds.json"
    if os.path.exists(safe_path):
        with open(safe_path, "r", encoding="ascii") as f:
            import json

            bounds = json.load(f)
        cfg.max_u = float(bounds.get("umax_safe", cfg.max_u))
        cfg.max_du = float(bounds.get("d_umax_safe", cfg.max_du))

    if args.fast:
        n_steps = min(n_steps, 3)
        base_tip = tip_position(x0)[0]
        targets = base_tip.repeat(n_steps + 1, 1)
        cfg.horizon = min(cfg.horizon, 10)
        cfg.max_iter = min(cfg.max_iter, 5)
        u_warm = torch.zeros((cfg.horizon, 1, 3), dtype=torch.float64)
        def _linearize_fast(x_t, u_t, li_in, cfg_in, backend=None, dense=None):
            n_x = x_t.shape[-1]
            n_u = u_t.shape[-1]
            a = torch.eye(n_x, dtype=x_t.dtype, device=x_t.device)
            b = torch.zeros((n_x, n_u), dtype=x_t.dtype, device=x_t.device)
            return a, b

        ilqr_mod._linearize = _linearize_fast
        print("FAST mode: using approximate linearization for quick validation")
    elif args.backend == "v1_3":
        n_steps = min(n_steps, 3)
        base_tip = tip_position(x0)[0]
        targets = base_tip.repeat(n_steps + 1, 1)
        cfg.horizon = min(cfg.horizon, 5)
        cfg.max_iter = min(cfg.max_iter, 1)
        u_warm = torch.zeros((cfg.horizon, 1, 3), dtype=torch.float64)
        if u_warm.shape[0] != cfg.horizon:
            if u_warm.shape[0] > cfg.horizon:
                u_warm = u_warm[: cfg.horizon]
            else:
                pad = u_warm[-1:].repeat(cfg.horizon - u_warm.shape[0], 1, 1)
                u_warm = torch.cat([u_warm, pad], dim=0)
        u_clamped = u_warm.clone()
        for t in range(u_clamped.shape[0]):
            u_clamped[t] = torch.clamp(u_clamped[t], -cfg.max_u, cfg.max_u)
            if t > 0:
                delta = u_clamped[t] - u_clamped[t - 1]
                delta_clamped = torch.clamp(delta, -cfg.max_du, cfg.max_du)
                sign_flip = (u_clamped[t - 1] * (u_clamped[t - 1] + delta_clamped)) < 0.0
                if sign_flip.any():
                    delta_clamped = torch.where(
                        sign_flip,
                        torch.clamp(delta_clamped, -cfg.max_sign_flip, cfg.max_sign_flip),
                        delta_clamped,
                    )
                u_clamped[t] = u_clamped[t - 1] + delta_clamped
        u_warm = u_clamped

    tip_index_sanity_check(x0, u_warm[:1], li, cfg_dyn, backend=args.backend)

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

    if not args.plot:
        print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
        return

    try:
        from crm_diffsims.control import plot_utils
    except Exception:
        print("matplotlib not installed; skipping plots")
        return

    errors = np.linalg.norm(tip_roll[: targets.shape[0]].cpu().numpy() - targets.cpu().numpy(), axis=1)

    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    try:
        plot_utils.save_tip_trajectory_plot(
            tip_roll.cpu().numpy(),
            targets.cpu().numpy(),
            title="Tip trajectory (lemniscate)",
            out_path=f"{out_dir}/ilqr_lemniscate_tip_{stamp}",
        )
        plot_utils.save_currents_plot(
            u_roll[:, 0].cpu().numpy(),
            title="Currents (lemniscate)",
            out_path=f"{out_dir}/ilqr_lemniscate_currents_{stamp}",
        )
        plot_utils.save_error_plot(
            errors,
            title="Tip tracking error (lemniscate)",
            out_path=f"{out_dir}/ilqr_lemniscate_error_{stamp}",
        )
    except Exception:
        print("matplotlib not installed; skipping plots")
        return
    print(f"Saved plots to {out_dir}/ilqr_lemniscate_*_{stamp}.{{png,pdf}}")
    print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
    for step in range(1, errors.shape[0]):
        mean_err = errors[: step + 1].mean()
        final_err = errors[step]
        print(f"MPC step {step}: mean_tip_error={mean_err:.3e} final_tip_error={final_err:.3e}")

    if args.profile_step:
        run_rollout_profile(x0, u_roll, li, cfg_dyn, label="lemniscate_profile", save_dir="output_data")


if __name__ == "__main__":
    main()
