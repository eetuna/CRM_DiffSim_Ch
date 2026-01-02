from datetime import datetime
import argparse
import os
import sys
import time

if __package__ is None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    python_root = os.path.join(repo_root, "python")
    if python_root not in sys.path:
        sys.path.insert(0, python_root)

import numpy as np
import torch

from crm_diffsims.control import ilqr as ilqr_mod
from crm_diffsims.control.ilqr import ILQRConfig, tip_index_sanity_check, tip_position
from crm_diffsims.control.run_report import write_run_report
from crm_diffsims.control.step_profiler import run_rollout_profile
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


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


def _atomic_save_npz(path, **payload):
    tmp_path = f"{path}.tmp.npz"
    np.savez(tmp_path, **payload)
    os.replace(tmp_path, path)


def _save_artifacts(
    out_path,
    targets,
    tip_xyz,
    u_roll,
    stats,
    li_mm,
    dt,
    umax_safe,
    d_umax_safe,
    dataset_u_used=None,
    dataset_currents_selected=None,
    dataset_targets_selected=None,
    n_total=None,
    start_idx=0,
    end_idx=-1,
    executed_tip_xyz=None,
    executed_u_roll=None,
    planned_first_tip_per_step=None,
    planned_terminal_tip_per_step=None,
    ramp_len=0,
    ramp_tip_xyz=None,
    ramp_u_applied=None,
    x0_after_ramp=None,
):
    tip_np = tip_xyz.detach().cpu().numpy()
    u_np = u_roll.detach().cpu().numpy() if u_roll is not None else np.zeros((0, 1, 3))
    target_np = targets.detach().cpu().numpy()
    errors = np.linalg.norm(tip_np[: target_np.shape[0]] - target_np[: tip_np.shape[0]], axis=1)
    costs_per_step = np.array([s.get("cost", np.nan) for s in stats], dtype=np.float64)
    terminal_err_per_step = errors[1 : len(stats) + 1]
    _atomic_save_npz(
        out_path,
        target_xyz=target_np,
        tip_xyz=tip_np,
        u_roll=u_np,
        dataset_u_used=dataset_u_used,
        executed_tip_xyz=executed_tip_xyz if executed_tip_xyz is not None else tip_np,
        executed_u_roll=executed_u_roll if executed_u_roll is not None else u_np,
        planned_first_tip_per_step=planned_first_tip_per_step,
        planned_terminal_tip_per_step=planned_terminal_tip_per_step,
        costs_per_step=costs_per_step,
        terminal_err_per_step=terminal_err_per_step,
        stats=np.array(stats, dtype=object),
        Li_mm=li_mm,
        dt=dt,
        umax_safe=umax_safe,
        d_umax_safe=d_umax_safe,
        dataset_currents_selected=dataset_currents_selected,
        dataset_targets_selected=dataset_targets_selected,
        N_total=n_total if n_total is not None else target_np.shape[0],
        start_idx=start_idx,
        end_idx=end_idx,
        ramp_len=ramp_len,
        ramp_tip_xyz=ramp_tip_xyz,
        ramp_u_applied=ramp_u_applied,
        x0_after_ramp=x0_after_ramp,
    )
    return errors


def _maybe_plot(out_dir, stamp, tip_np, target_np, u_np, errors, planned_np=None):
    try:
        from crm_diffsims.control import plot_utils
    except Exception:
        print("matplotlib not installed; skipping plots")
        return
    plot_utils.save_tip_trajectory_plot(
        tip_np,
        target_np,
        title="Tip trajectory (circle)",
        out_path=f"{out_dir}/ilqr_circle_tip_{stamp}",
        planned_xyz=planned_np,
    )
    if u_np.shape[0]:
        plot_utils.save_currents_plot(
            u_np[:, 0],
            title="Currents (circle)",
            out_path=f"{out_dir}/ilqr_circle_currents_{stamp}",
        )
    plot_utils.save_error_plot(
        errors,
        title="Tip tracking error (circle)",
        out_path=f"{out_dir}/ilqr_circle_error_{stamp}",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast", action="store_true", help="Run a short validation rollout")
    parser.add_argument("--backend", choices=["v1_3", "autograd"], default="v1_3")
    parser.add_argument("--plot", action="store_true", help="Save plots to output_data/")
    parser.add_argument("--profile-step", action="store_true", help="Profile crm_step and save flight recorder")
    parser.add_argument("--max-wall-sec", type=float, default=300.0, help="Hard wall time limit in seconds")
    parser.add_argument("--save-every", type=int, default=1, help="Save NPZ every N MPC steps")
    parser.add_argument("--heartbeat-sec", type=float, default=5.0, help="Print progress every N seconds")
    parser.add_argument("--init-u", choices=["zero", "dataset", "ramp_to_dataset"], default="ramp_to_dataset")
    parser.add_argument("--bounds-mode", choices=["fixed", "ramp"], default="fixed")
    parser.add_argument("--umax-start", type=float, default=0.1)
    parser.add_argument("--umax-max", type=float, default=0.3)
    parser.add_argument("--dumax-start", type=float, default=0.01)
    parser.add_argument("--dumax-max", type=float, default=0.05)
    parser.add_argument("--ramp-steps", type=int, default=30)
    parser.add_argument("--profile-lite", action="store_true")
    parser.add_argument("--replay-dataset-mode", action="store_true")
    parser.add_argument("--replay-only", action="store_true")
    parser.add_argument("--mode", choices=["replay_only_full", "replay_dataset_mode_full"], default=None)
    parser.add_argument("--mpc-steps", type=int, default=None)
    parser.add_argument("--horizon", type=int, default=None)
    parser.add_argument("--max-iter", type=int, default=None)
    parser.add_argument("--start-idx", type=int, default=0)
    parser.add_argument("--end-idx", type=int, default=-1)
    parser.add_argument("--full-dataset", action="store_true")
    parser.add_argument("--init-from-ramp", action="store_true")
    parser.add_argument("--ramp-len", type=int, default=120)
    parser.add_argument("--plot-every", type=int, default=10)
    parser.add_argument("--save-u-seq", action="store_true", help="Save iLQR u_seq after MPC run")
    args = parser.parse_args()

    if args.mode == "replay_only_full":
        args.replay_only = True
        args.full_dataset = True
    elif args.mode == "replay_dataset_mode_full":
        args.replay_dataset_mode = True
        args.full_dataset = True

    if args.replay_dataset_mode or args.replay_only:
        args.init_u = "dataset"
        args.bounds_mode = "fixed"
        if args.mpc_steps is None:
            args.mpc_steps = 30
        if args.horizon is None:
            args.horizon = 15
        if args.max_iter is None:
            args.max_iter = 3

    if args.profile_lite:
        args.bounds_mode = "ramp"
        args.max_wall_sec = 120.0

    torch.manual_seed(0)
    np.random.seed(0)
    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)

    dataset_path = "data/dyn_fk_ramp_circle1_hold1.npz"
    target_np, currents_np, dt, li_val = _load_targets(dataset_path)

    n_steps = 80
    targets = torch.tensor(target_np[: n_steps + 1], dtype=torch.float64)
    x0, _ = build_state()
    li = torch.tensor([li_val], dtype=torch.float64)

    horizon = 30
    if args.horizon is not None:
        horizon = args.horizon

    cfg = ILQRConfig(
        horizon=horizon,
        max_iter=10,
        w_tip=1.0,
        w_tip_terminal=25.0,
        w_u=5e-3,
        w_du=5e-2,
        max_du=0.008,
        max_sign_flip=0.004,
    )
    cfg.linearization_backend = args.backend
    cfg.linearization_dense = True

    if args.profile_lite:
        n_steps = min(n_steps, 10)
        targets = targets[: n_steps + 1]
        cfg.max_iter = 5
    if args.full_dataset:
        full_targets_np = target_np
        full_currents_np = currents_np if currents_np is not None else np.zeros((0, 3))
        end_idx = full_targets_np.shape[0] if args.end_idx < 0 else min(args.end_idx, full_targets_np.shape[0])
        start_idx = max(0, min(args.start_idx, end_idx))
        args.start_idx = start_idx
        args.end_idx = end_idx
        target_np = full_targets_np[start_idx:end_idx]
        currents_np = full_currents_np[start_idx:end_idx] if full_currents_np.shape[0] else full_currents_np
        targets = torch.tensor(target_np, dtype=torch.float64)
        n_steps = max(0, targets.shape[0] - 1)
        args.mpc_steps = n_steps

    dataset_targets_selected_full = np.asarray(target_np, dtype=np.float64)
    if currents_np is not None and currents_np.shape[0]:
        dataset_currents_selected_full = np.asarray(currents_np, dtype=np.float64)
    else:
        dataset_currents_selected_full = np.zeros((dataset_targets_selected_full.shape[0], 3), dtype=np.float64)

    ramp_len = 0
    ramp_tip_xyz = None
    ramp_u_applied = None
    x0_after_ramp = None
    if args.init_from_ramp and currents_np is not None and currents_np.shape[0]:
        ramp_len = min(args.ramp_len, currents_np.shape[0], target_np.shape[0])
        if ramp_len > 0:
            u_ramp = torch.tensor(currents_np[:ramp_len], dtype=torch.float64).view(ramp_len, 1, 3)
            x_ramp = x0.clone()
            ramp_tip = [tip_position(x_ramp)[0].detach().cpu().numpy()]
            ramp_u = []
            for t in range(ramp_len):
                u_t = u_ramp[t : t + 1]
                ramp_u.append(u_t.detach().cpu().numpy())
                x_ramp, _ = crm_step_v1_3_forward(x_ramp, u_t, li, cfg_dyn)
                ramp_tip.append(tip_position(x_ramp)[0].detach().cpu().numpy())
            x0 = x_ramp
            x0_after_ramp = x_ramp.detach().cpu().numpy()
            ramp_tip_xyz = np.asarray(ramp_tip, dtype=np.float64)
            ramp_u_applied = np.asarray(ramp_u, dtype=np.float64)
            target_np = target_np[ramp_len:]
            if currents_np is not None and currents_np.shape[0]:
                currents_np = currents_np[ramp_len:]
            targets = torch.tensor(target_np, dtype=torch.float64)
            n_steps = max(0, targets.shape[0] - 1)

    if currents_np is not None:
        u_seed_full = torch.tensor(currents_np, dtype=torch.float64).view(-1, 1, 3)
        u_seed = u_seed_full[:horizon]
    else:
        u_seed_full = None
        u_seed = None
    u_warm = torch.zeros((horizon, 1, 3), dtype=torch.float64)
    if args.mpc_steps is not None:
        n_steps = min(n_steps, args.mpc_steps)
        targets = targets[: n_steps + 1]
    if args.max_iter is not None:
        cfg.max_iter = args.max_iter

    if args.replay_dataset_mode or args.replay_only:
        cfg.max_u = 0.35
        cfg.max_du = 0.05
        cfg.use_safe_bounds = False
        safe_max_u = cfg.max_u
        safe_max_du = cfg.max_du
        if args.replay_only:
            print(f"REPLAY_DATASET_MODE=1 init_u=dataset umax=0.35 d_umax=0.05 Li_mm={li_val:.1f}")
            print("REPLAY_ONLY=1")
        else:
            print(f"REPLAY_DATASET_MODE=1 init_u=dataset umax=0.35 d_umax=0.05 Li_mm={li_val:.1f}")
    else:
        safe_path = "docs/control/safe_bounds.json"
        if os.path.exists(safe_path):
            with open(safe_path, "r", encoding="ascii") as f:
                import json

                bounds = json.load(f)
            cfg.max_u = float(bounds.get("umax_safe", cfg.max_u))
            cfg.max_du = float(bounds.get("d_umax_safe", cfg.max_du))
        safe_max_u = cfg.max_u
        safe_max_du = cfg.max_du
        if args.bounds_mode == "ramp":
            cfg.max_u = min(safe_max_u, args.umax_start)
            cfg.max_du = min(safe_max_du, args.dumax_start)
        print(
            f"init_u={args.init_u}, bounds_mode={args.bounds_mode}, Li_mm={li_val:.1f}, "
            f"umax_safe={safe_max_u:.4f}, d_umax_safe={safe_max_du:.4f}"
        )

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
    if u_warm.shape[0] != cfg.horizon:
        if u_warm.shape[0] > cfg.horizon:
            u_warm = u_warm[: cfg.horizon]
        else:
            pad = u_warm[-1:].repeat(cfg.horizon - u_warm.shape[0], 1, 1)
            u_warm = torch.cat([u_warm, pad], dim=0)
    if args.init_u == "dataset" and u_seed is not None:
        u_warm = u_seed.clone()
    elif args.init_u == "ramp_to_dataset" and u_seed is not None:
        u0 = u_warm[0:1].clone()
        ramped = []
        n_ramp = max(1, u_seed.shape[0])
        for t in range(u_seed.shape[0]):
            alpha = min(1.0, float(t + 1) / n_ramp)
            ramped.append((1.0 - alpha) * u0 + alpha * u_seed[t : t + 1])
        u_warm = torch.cat(ramped, dim=0)
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

    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.replay_dataset_mode or args.replay_only:
        out_path = f"{out_dir}/replay_mode_circle_{stamp}.npz"
    else:
        out_path = f"{out_dir}/ilqr_circle_{stamp}.npz"
    start_time = time.time()
    last_heartbeat = start_time
    last_result = None
    max_wall_hit = False

    if not args.replay_dataset_mode and not args.replay_only:
        ilqr_mod._maybe_apply_safe_bounds(cfg)
    n_mpc_steps = targets.shape[0] - 1
    u_seq = u_warm.clone()
    x_roll = [x0]
    u_roll = []
    stats = []
    dataset_u_used = []
    planned_first_tip_per_step = []
    planned_terminal_tip_per_step = []
    planned_last_tip = None
    dataset_targets_selected = dataset_targets_selected_full
    n_total = dataset_targets_selected.shape[0]
    start_idx = args.start_idx
    end_idx = args.end_idx
    dataset_currents_selected = dataset_currents_selected_full

    def _emit_report(report_errors):
        mean_err = float(np.mean(report_errors)) if report_errors.size else float("nan")
        max_err = float(np.max(report_errors)) if report_errors.size else float("nan")
        report_path = os.path.join(
            "docs/control/run_reports", f"run_ilqr_circle_{stamp}.md"
        )
        write_run_report(
            report_path=report_path,
            command=" ".join(sys.argv),
            dataset_name=os.path.basename(dataset_path),
            start_idx=start_idx,
            end_idx=end_idx,
            n_total=n_total,
            li_mm=li_val,
            dt=dt,
            umax=cfg.max_u,
            d_umax=cfg.max_du,
            mean_error=mean_err,
            max_error=max_err,
            artifact_npz=out_path,
            plot_prefix=f"{out_dir}/ilqr_circle_{stamp}" if args.plot else None,
            init_from_ramp=args.init_from_ramp,
            ramp_len=ramp_len,
            max_wall_hit=max_wall_hit,
        )
    tip_roll = tip_position(torch.cat(x_roll, dim=0))
    errors = _save_artifacts(
        out_path,
        targets,
        tip_roll,
        None,
        stats,
        li_val,
        dt,
        cfg.max_u,
        cfg.max_du,
        dataset_u_used=np.asarray(dataset_u_used, dtype=np.float64).reshape(-1, 1, 3) if dataset_u_used else None,
        dataset_currents_selected=dataset_currents_selected,
        dataset_targets_selected=dataset_targets_selected,
        n_total=n_total,
        start_idx=start_idx,
        end_idx=end_idx,
        executed_tip_xyz=tip_roll.detach().cpu().numpy(),
        executed_u_roll=None,
        planned_first_tip_per_step=np.asarray(planned_first_tip_per_step, dtype=np.float64),
        planned_terminal_tip_per_step=np.asarray(planned_terminal_tip_per_step, dtype=np.float64),
        ramp_len=ramp_len,
        ramp_tip_xyz=ramp_tip_xyz,
        ramp_u_applied=ramp_u_applied,
        x0_after_ramp=x0_after_ramp,
    )
    if args.plot:
        _maybe_plot(
            out_dir,
            stamp,
            tip_roll.detach().cpu().numpy(),
            targets.detach().cpu().numpy(),
            np.zeros((0, 1, 3)),
            errors,
            planned_np=planned_last_tip,
        )

    for t_step in range(n_mpc_steps):
        if args.bounds_mode == "ramp":
            alpha = min(1.0, max(0.0, t_step / max(1, args.ramp_steps)))
            cfg.max_u = min(safe_max_u, args.umax_start + (args.umax_max - args.umax_start) * alpha)
            cfg.max_du = min(safe_max_du, args.dumax_start + (args.dumax_max - args.dumax_start) * alpha)
        if (args.replay_dataset_mode or args.replay_only) and u_seed_full is not None:
            u_window = u_seed_full[t_step : t_step + cfg.horizon]
            if u_window.shape[0] < cfg.horizon:
                pad = u_window[-1:].repeat(cfg.horizon - u_window.shape[0], 1, 1)
                u_window = torch.cat([u_window, pad], dim=0)
            u_seq = u_window.clone()
            dataset_u_used.append(u_window[0].detach().cpu().numpy())
        horizon_targets = targets[t_step : t_step + cfg.horizon + 1]
        if horizon_targets.shape[0] < cfg.horizon + 1:
            pad = horizon_targets[-1:].repeat(cfg.horizon + 1 - horizon_targets.shape[0], 1)
            horizon_targets = torch.cat([horizon_targets, pad], dim=0)
        if args.replay_only:
            result = None
            stats.append(
                {
                    "clamp_hits": 0,
                    "sign_flip_hits": 0,
                    "clamp_rate": 0.0,
                    "sign_flip_rate": 0.0,
                    "nan_detected": False,
                    "first_iter_u_seq": None,
                    "timed_out": False,
                    "cost": float("nan"),
                    "converged": True,
                }
            )
        else:
            if args.replay_dataset_mode and args.init_u == "dataset" and u_seed_full is not None:
                u_seq_head = u_seq[:3].detach().cpu().numpy()
                u_ds_head = u_seed_full[t_step : t_step + 3].detach().cpu().numpy()
                print(f"u_seq_head={u_seq_head} u_dataset_head={u_ds_head}")
                if not np.allclose(u_seq_head, u_ds_head, atol=1e-12):
                    raise RuntimeError("Replay dataset warm-start mismatch: u_seq != dataset window")
            result = ilqr_mod.ilqr_solve(
                x_roll[-1],
                u_seq,
                horizon_targets,
                li,
                cfg_dyn,
                cfg,
                max_wall_sec=args.max_wall_sec,
                start_time=start_time,
            )
            last_result = result
            if result.x_seq is not None:
                planned_tip = tip_position(result.x_seq).detach().cpu().numpy()
                planned_first_tip_per_step.append(planned_tip[0])
                planned_terminal_tip_per_step.append(planned_tip[-1])
                planned_last_tip = planned_tip
            stats.append(
                {
                    "clamp_hits": result.clamp_hits,
                    "sign_flip_hits": result.sign_flip_hits,
                    "clamp_rate": result.clamp_rate,
                    "sign_flip_rate": result.sign_flip_rate,
                    "nan_detected": result.nan_detected,
                    "first_iter_u_seq": result.first_iter_u_seq,
                    "timed_out": result.timed_out,
                    "cost": result.cost,
                    "converged": result.converged,
                }
            )
        print(
            f"MPC step {t_step + 1}/{n_mpc_steps}: "
            f"umax={cfg.max_u:.4f} dumax={cfg.max_du:.4f} "
            f"clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}"
        )
        if args.replay_only:
            u_apply = u_seq[0:1]
        else:
            u_apply = result.u_seq[0:1]
            if args.replay_dataset_mode and result.nan_detected and u_seed_full is not None:
                u_apply = u_seq[0:1]
        u_roll.append(u_apply)
        if (args.replay_dataset_mode or args.replay_only) and u_seed_full is not None:
            u_dataset0 = u_seq[0:1].detach().cpu().numpy()
            u_applied = u_apply.detach().cpu().numpy()
            u_diff = np.linalg.norm(u_applied - u_dataset0)
            print(
                "replay_step: "
                f"u_dataset0={u_dataset0} u_applied={u_applied} "
                f"diff_norm={u_diff:.3e} umax={cfg.max_u:.4f} dumax={cfg.max_du:.4f}"
            )
        x_next, step_ok = ilqr_mod._step_forward(
            x_roll[-1],
            u_apply,
            li,
            cfg_dyn,
            cfg.linearization_backend,
        )
        if not step_ok or not torch.isfinite(x_next).all().item():
            break
        x_roll.append(x_next)
        if not args.replay_only:
            u_seq = torch.cat([result.u_seq[1:], result.u_seq[-1:]], dim=0)
        tip_roll = tip_position(torch.cat(x_roll, dim=0))
        u_roll_cat = torch.cat(u_roll, dim=0) if u_roll else None
        if t_step % max(1, args.save_every) == 0 or t_step == n_mpc_steps - 1:
            errors = _save_artifacts(
                out_path,
                targets,
                tip_roll,
                u_roll_cat,
                stats,
                li_val,
                dt,
                cfg.max_u,
                cfg.max_du,
                dataset_u_used=np.asarray(dataset_u_used, dtype=np.float64).reshape(-1, 1, 3) if dataset_u_used else None,
                dataset_currents_selected=dataset_currents_selected,
                dataset_targets_selected=dataset_targets_selected,
                n_total=n_total,
                start_idx=start_idx,
                end_idx=end_idx,
                executed_tip_xyz=tip_roll.detach().cpu().numpy(),
                executed_u_roll=u_roll_cat.detach().cpu().numpy() if u_roll_cat is not None else None,
                planned_first_tip_per_step=np.asarray(planned_first_tip_per_step, dtype=np.float64),
                planned_terminal_tip_per_step=np.asarray(planned_terminal_tip_per_step, dtype=np.float64),
                ramp_len=ramp_len,
                ramp_tip_xyz=ramp_tip_xyz,
                ramp_u_applied=ramp_u_applied,
                x0_after_ramp=x0_after_ramp,
            )
            if args.plot and (t_step % max(1, args.plot_every) == 0 or t_step == n_mpc_steps - 1):
                _maybe_plot(
                    out_dir,
                    stamp,
                    tip_roll.detach().cpu().numpy(),
                    targets.detach().cpu().numpy(),
                    u_roll_cat.detach().cpu().numpy() if u_roll_cat is not None else np.zeros((0, 1, 3)),
                    errors,
                    planned_np=planned_last_tip,
                )
        now = time.time()
        if now - last_heartbeat >= args.heartbeat_sec:
            last_heartbeat = now
            elapsed = now - start_time
            last_cost = stats[-1]["cost"] if stats else float("nan")
            last_err = errors[-1] if errors.size else float("nan")
            last_clamp = stats[-1]["clamp_rate"] if stats else float("nan")
            last_sign = stats[-1]["sign_flip_rate"] if stats else float("nan")
            iter_count = last_result.iterations if last_result is not None else 0
            print(
                "heartbeat: "
                f"elapsed={elapsed:.1f}s "
                f"step={t_step + 1}/{n_mpc_steps} "
                f"ilqr_iter={iter_count} "
                f"cost={last_cost:.3e} "
                f"terminal_err={last_err:.3e} "
                f"clamp_rate={last_clamp:.2%} "
                f"sign_flip_rate={last_sign:.2%}"
            )
        if args.max_wall_sec is not None and (now - start_time) > args.max_wall_sec:
            print(f"Max wall time exceeded ({args.max_wall_sec:.1f}s); saving partial results.")
            max_wall_hit = True
            errors = _save_artifacts(
                out_path,
                targets,
                tip_roll,
                u_roll_cat,
                stats,
                li_val,
                dt,
                cfg.max_u,
                cfg.max_du,
                dataset_u_used=np.asarray(dataset_u_used, dtype=np.float64).reshape(-1, 1, 3) if dataset_u_used else None,
                dataset_currents_selected=dataset_currents_selected,
                dataset_targets_selected=dataset_targets_selected,
                n_total=n_total,
                start_idx=start_idx,
                end_idx=end_idx,
                executed_tip_xyz=tip_roll.detach().cpu().numpy(),
                executed_u_roll=u_roll_cat.detach().cpu().numpy() if u_roll_cat is not None else None,
                planned_first_tip_per_step=np.asarray(planned_first_tip_per_step, dtype=np.float64),
                planned_terminal_tip_per_step=np.asarray(planned_terminal_tip_per_step, dtype=np.float64),
                ramp_len=ramp_len,
                ramp_tip_xyz=ramp_tip_xyz,
                ramp_u_applied=ramp_u_applied,
                x0_after_ramp=x0_after_ramp,
            )
            _emit_report(errors)
            if args.plot:
                _maybe_plot(
                    out_dir,
                    stamp,
                    tip_roll.detach().cpu().numpy(),
                    targets.detach().cpu().numpy(),
                    u_roll_cat.detach().cpu().numpy() if u_roll_cat is not None else np.zeros((0, 1, 3)),
                    errors,
                    planned_np=planned_last_tip,
                )
            sys.exit(2)
        if not args.replay_only:
            if result.timed_out:
                break
            if result.nan_detected and not args.replay_dataset_mode:
                break

    tip_roll = tip_position(torch.cat(x_roll, dim=0))
    u_roll_cat = torch.cat(u_roll, dim=0) if u_roll else None
    errors = _save_artifacts(
        out_path,
        targets,
        tip_roll,
        u_roll_cat,
        stats,
        li_val,
        dt,
        cfg.max_u,
        cfg.max_du,
        dataset_u_used=np.asarray(dataset_u_used, dtype=np.float64).reshape(-1, 1, 3) if dataset_u_used else None,
        dataset_currents_selected=dataset_currents_selected,
        dataset_targets_selected=dataset_targets_selected,
        n_total=n_total,
        start_idx=start_idx,
        end_idx=end_idx,
        executed_tip_xyz=tip_roll.detach().cpu().numpy(),
        executed_u_roll=u_roll_cat.detach().cpu().numpy() if u_roll_cat is not None else None,
        planned_first_tip_per_step=np.asarray(planned_first_tip_per_step, dtype=np.float64),
        planned_terminal_tip_per_step=np.asarray(planned_terminal_tip_per_step, dtype=np.float64),
        ramp_len=ramp_len,
        ramp_tip_xyz=ramp_tip_xyz,
        ramp_u_applied=ramp_u_applied,
        x0_after_ramp=x0_after_ramp,
    )
    _emit_report(errors)
    if args.plot:
        _maybe_plot(
            out_dir,
            stamp,
            tip_roll.detach().cpu().numpy(),
            targets.detach().cpu().numpy(),
            u_roll_cat.detach().cpu().numpy() if u_roll_cat is not None else np.zeros((0, 1, 3)),
            errors,
            planned_np=planned_last_tip,
        )

    if args.save_u_seq:
        first_iter_u = None
        if stats and stats[0].get("first_iter_u_seq") is not None:
            first_iter_u = stats[0]["first_iter_u_seq"].detach().cpu().numpy()
        np.savez(
            f"{out_dir}/ilqr_circle_u_seq_{stamp}.npz",
            u_seq=u_seq.detach().cpu().numpy(),
            first_iter_u_seq=first_iter_u,
        )

    if any(s.get("nan_detected") for s in stats) and not args.replay_dataset_mode:
        _emit_report(errors)
        print("NaN detected in iLQR; saved diagnostics.")
        sys.exit(1)

    if args.fast:
        print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
        return

    if not args.plot:
        print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
        return

    print(f"Saved plots to {out_dir}/ilqr_circle_*_{stamp}.{{png,pdf}}")
    print(f"Final MPC clamp_rate={stats[-1]['clamp_rate']:.2%} sign_flip_rate={stats[-1]['sign_flip_rate']:.2%}")
    for step in range(1, errors.shape[0]):
        mean_err = errors[: step + 1].mean()
        final_err = errors[step]
        print(f"MPC step {step}: mean_tip_error={mean_err:.3e} final_tip_error={final_err:.3e}")

    if args.profile_step:
        if u_roll_cat is not None and u_roll_cat.shape[0] > 0:
            run_rollout_profile(x0, u_roll_cat, li, cfg_dyn, label="circle_profile", save_dir="output_data")


if __name__ == "__main__":
    main()
