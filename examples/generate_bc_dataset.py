import argparse
import io
import json
import os
import sys
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime

if __package__ is None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    python_root = os.path.join(repo_root, "python")
    if python_root not in sys.path:
        sys.path.insert(0, python_root)

import numpy as np
import torch

from crm_diffsims.control.cem_mpc import CEMConfig, cem_mpc
from crm_diffsims.control.ilqr import ILQRConfig, run_mpc
from crm_diffsims.control.run_report import write_run_report
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


TIP_IDX = slice(24, 27)


def _load_safe_bounds():
    path = os.path.join("docs", "control", "safe_bounds.json")
    if not os.path.exists(path):
        return 0.1, 0.01
    with open(path, "r", encoding="ascii") as f:
        payload = json.load(f)
    return float(payload.get("umax_safe", 0.1)), float(payload.get("d_umax_safe", 0.01))


def _load_dataset(path):
    data = np.load(path)
    if "tip_desired" in data:
        candidate = data["tip_desired"]
        if np.isfinite(candidate).all():
            target = candidate
        elif "tip_fk" in data:
            target = data["tip_fk"]
        else:
            target = data["tip_dyn"]
    elif "tip_fk" in data:
        target = data["tip_fk"]
    else:
        target = data["tip_dyn"]
    currents = data.get("currents")
    dt = float(data.get("dt", 0.2))
    li = float(data.get("insertion_length", 94.3))
    segments = data.get("segments")
    return target, currents, dt, li, segments


def _infer_ramp_len(segments, args_ramp_len):
    if args_ramp_len is not None:
        return int(args_ramp_len)
    if segments is None:
        return 0
    if np.ndim(segments) == 0:
        return int(segments)
    if len(segments) > 0:
        return int(segments[0])
    return 0


def _rate_limit(u_seq, d_umax):
    u_limited = u_seq.copy()
    for t in range(1, u_limited.shape[0]):
        delta = u_limited[t] - u_limited[t - 1]
        delta = np.clip(delta, -d_umax, d_umax)
        u_limited[t] = u_limited[t - 1] + delta
    return u_limited


def _clamp_currents(u_seq, umax, d_umax):
    u = np.clip(u_seq, -umax, umax)
    return _rate_limit(u, d_umax)


def _rollout_prefix(x0, u_prefix, li, cfg_dyn):
    if u_prefix.shape[0] == 0:
        return x0
    x_t = x0
    for t in range(u_prefix.shape[0]):
        u_t = torch.tensor(u_prefix[t : t + 1], dtype=torch.float64).view(1, 1, 3)
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        solver_exit = int(cache["solver_exit"].item())
        residual = float(cache["residual_norm"].item())
        if solver_exit != 0 or residual > cfg_dyn.residual_threshold:
            raise RuntimeError(f"nonconvergence during prefix at step {t}")
        if not torch.isfinite(x_tp1).all().item():
            raise RuntimeError(f"nonfinite during prefix at step {t}")
        x_t = x_tp1
    return x_t


def _rollout_collect(x0, u_seq, targets, li, cfg_dyn):
    x_list = []
    tip_list = []
    solver_exit = []
    residual_norm = []
    unbounded = []
    nonfinite = []

    x_t = x0
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
        for t in range(u_seq.shape[0]):
            x_list.append(x_t.detach().cpu().numpy()[0])
            tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())
            u_t = u_seq[t : t + 1]
            x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
            exit_flag = int(cache["solver_exit"].item())
            residual = float(cache["residual_norm"].item())
            solver_exit.append(exit_flag)
            residual_norm.append(residual)
            message = stdout_buf.getvalue() + stderr_buf.getvalue()
            unbounded_flag = "Unbounded" in message
            unbounded.append(unbounded_flag)
            nonfinite_flag = not torch.isfinite(x_tp1).all().item()
            nonfinite.append(nonfinite_flag)

            if exit_flag != 0 or residual > cfg_dyn.residual_threshold:
                break
            if unbounded_flag or nonfinite_flag:
                break
            x_t = x_tp1

    tip_np = np.asarray(tip_list, dtype=np.float64)
    target_np = targets[: tip_np.shape[0]]
    errors = np.linalg.norm(tip_np - target_np, axis=1) if tip_np.size else np.zeros((0,))

    return (
        np.asarray(x_list, dtype=np.float64),
        tip_np,
        np.asarray(solver_exit, dtype=np.int32),
        np.asarray(residual_norm, dtype=np.float64),
        np.asarray(unbounded, dtype=bool),
        np.asarray(nonfinite, dtype=bool),
        errors,
    )


def _prepare_targets(targets, start_idx, n_steps, wrap_len=None):
    if wrap_len:
        idx = (start_idx + np.arange(n_steps + 1)) % wrap_len
        return targets[idx]
    end_idx = start_idx + n_steps + 1
    clipped = targets[start_idx:end_idx]
    if clipped.shape[0] < n_steps + 1:
        pad = np.repeat(clipped[-1:], n_steps + 1 - clipped.shape[0], axis=0)
        clipped = np.concatenate([clipped, pad], axis=0)
    return clipped


def _slice_with_wrap(arr, start_idx, n_steps, wrap_len=None):
    if wrap_len:
        idx = (start_idx + np.arange(n_steps)) % wrap_len
        return arr[idx]
    end_idx = start_idx + n_steps
    clipped = arr[start_idx:end_idx]
    if clipped.shape[0] < n_steps:
        pad = np.repeat(clipped[-1:], n_steps - clipped.shape[0], axis=0)
        clipped = np.concatenate([clipped, pad], axis=0)
    return clipped


def _episode_indices(start_idx, end_idx, episode_len, episodes):
    if episodes is None:
        total = max(end_idx - start_idx, 0)
        episodes = max(int(np.ceil(total / episode_len)), 1)
    for ep in range(episodes):
        s = start_idx + ep * episode_len
        e = s + episode_len
        if end_idx is not None and s >= end_idx:
            break
        yield ep, s, e


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, default="data/dyn_fk_ramp_circle1_hold1.npz")
    parser.add_argument(
        "--controller",
        choices=["replay_only_full", "cem", "cem_mpc", "ilqr_warmstart"],
        default="replay_only_full",
    )
    parser.add_argument("--preset", choices=["circle_easy", "lemniscate_staged", "mixed_easy"], default=None)
    parser.add_argument("--episodes", type=int, default=1)
    parser.add_argument("--episode-len", type=int, default=60)
    parser.add_argument("--start-idx", type=int, default=0)
    parser.add_argument("--end-idx", type=int, default=-1)
    parser.add_argument("--init-from-ramp", action="store_true")
    parser.add_argument("--ramp-len", type=int, default=None)
    parser.add_argument("--wrap-dataset", action="store_true")
    parser.add_argument(
        "--filter-policy",
        choices=["keep_all", "drop_unbounded", "drop_any_fail"],
        default="drop_unbounded",
    )
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-npz", type=str, default=None)
    parser.add_argument("--report-path", type=str, default=None)
    parser.add_argument("--max-wall-sec", type=float, default=120.0)
    parser.add_argument("--backend", choices=["v1_3", "autograd"], default="v1_3")
    parser.add_argument("--horizon", type=int, default=20)
    parser.add_argument("--mpc-steps", type=int, default=None)
    parser.add_argument("--max-iter", type=int, default=5)
    parser.add_argument("--num-samples", type=int, default=32)
    parser.add_argument("--num-elites", type=int, default=6)
    parser.add_argument("--cem-iters", type=int, default=2)
    parser.add_argument("--u-std-init", type=float, default=0.02)
    parser.add_argument("--umax", type=float, default=None)
    parser.add_argument("--d-umax", type=float, default=None)
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    preset_datasets = None
    if args.preset == "circle_easy":
        args.controller = "replay_only_full"
        args.init_from_ramp = True
        preset_datasets = [("circle", "data/dyn_fk_ramp_circle1_hold1.npz", args.episodes)]
    elif args.preset == "lemniscate_staged":
        args.controller = "replay_only_full"
        args.init_from_ramp = True
        preset_datasets = [("lemniscate", "data/dyn_fk_lem1_y40_a10_L94_hold1.npz", args.episodes)]
    elif args.preset == "mixed_easy":
        args.controller = "replay_only_full"
        args.init_from_ramp = True
        half = args.episodes // 2
        preset_datasets = [
            ("circle", "data/dyn_fk_ramp_circle1_hold1.npz", half),
            ("lemniscate", "data/dyn_fk_lem1_y40_a10_L94_hold1.npz", args.episodes - half),
        ]

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0_default, li_default = build_state()

    umax_safe, d_umax_safe = _load_safe_bounds()
    if args.umax is not None:
        umax_safe = float(args.umax)
    if args.d_umax is not None:
        d_umax_safe = float(args.d_umax)

    all_x = []
    all_target = []
    all_u = []
    all_episode = []
    all_step = []
    all_solver = []
    all_residual = []
    all_unbounded = []
    all_tip = []
    all_nonfinite = []
    all_x_raw = []
    all_target_raw = []
    all_u_raw = []
    all_episode_raw = []
    all_step_raw = []
    all_solver_raw = []
    all_residual_raw = []
    all_unbounded_raw = []
    all_tip_raw = []
    all_nonfinite_raw = []
    all_keep_raw = []
    all_drop_reason_raw = []
    dataset_names = []
    total_steps = 0
    kept_steps = 0
    dropped_steps = 0
    drop_reasons = {
        "unbounded": 0,
        "solver_exit": 0,
        "residual": 0,
        "nonfinite": 0,
        "target_nonfinite": 0,
        "policy_drop": 0,
        "prefix_fail": 0,
    }

    max_wall_hit = False
    total_n_total = 0
    li_report = 0.0
    dt_report = 0.0
    ramp_len_report = 0
    episode_offset = 0

    dataset_list = preset_datasets or [(None, args.dataset, args.episodes)]
    for dataset_label, dataset_path, episodes in dataset_list:
        if episodes <= 0:
            continue
        targets, currents, dt, li_val, segments = _load_dataset(dataset_path)
        li = torch.tensor([li_val], dtype=torch.float64)
        dataset_names.append(os.path.basename(dataset_path))
        total_n_total += targets.shape[0]
        if li_report == 0.0:
            li_report = li_val
        if dt_report == 0.0:
            dt_report = dt

        n_total = targets.shape[0]
        if args.wrap_dataset:
            end_idx = None
        else:
            end_idx = n_total if args.end_idx < 0 else int(args.end_idx)
        start_idx = max(0, args.start_idx)

        ramp_len = _infer_ramp_len(segments, args.ramp_len) if args.init_from_ramp else 0
        ramp_len_report = max(ramp_len_report, ramp_len)

        if args.controller == "replay_only_full" and currents is None:
            raise SystemExit("replay_only_full requires currents in the dataset")

        if currents is None:
            currents = np.zeros((n_total, 3), dtype=np.float64)

        for ep_id, s_idx, e_idx in _episode_indices(start_idx, end_idx, args.episode_len, episodes):
            n_steps = e_idx - s_idx
            if n_steps <= 0:
                continue

            x0 = x0_default.clone()
            try:
                s_idx_prefix = s_idx % n_total if args.wrap_dataset else s_idx
                if args.init_from_ramp and ramp_len > 0:
                    ramp_u = _clamp_currents(currents[:ramp_len], umax_safe, d_umax_safe)
                    x0 = _rollout_prefix(x0, ramp_u, li, cfg_dyn)
                    if s_idx_prefix > ramp_len:
                        prefix_raw = _slice_with_wrap(currents, ramp_len, s_idx_prefix - ramp_len)
                        prefix_u = _clamp_currents(prefix_raw, umax_safe, d_umax_safe)
                        x0 = _rollout_prefix(x0, prefix_u, li, cfg_dyn)
                elif s_idx_prefix > 0:
                    prefix_raw = _slice_with_wrap(currents, 0, s_idx_prefix)
                    prefix_u = _clamp_currents(prefix_raw, umax_safe, d_umax_safe)
                    x0 = _rollout_prefix(x0, prefix_u, li, cfg_dyn)
            except RuntimeError as exc:
                drop_reasons["prefix_fail"] += 1
                print(f"prefix failed for episode {ep_id}: {exc}")
                continue

            wrap_len = n_total if args.wrap_dataset else None
            target_seq = _prepare_targets(targets, s_idx, n_steps, wrap_len=wrap_len)

            if args.controller == "replay_only_full":
                u_raw = _slice_with_wrap(currents, s_idx, n_steps, wrap_len=wrap_len)
                u_seq_np = _clamp_currents(u_raw, umax_safe, d_umax_safe)
                u_seq = torch.tensor(u_seq_np, dtype=torch.float64).view(n_steps, 1, 3)
            elif args.controller in ("cem", "cem_mpc"):
                cfg = CEMConfig(
                    horizon=args.horizon,
                    mpc_steps=n_steps if args.mpc_steps is None else args.mpc_steps,
                    num_samples=args.num_samples,
                    num_elites=args.num_elites,
                    cem_iters=args.cem_iters,
                    u_std_init=args.u_std_init,
                    umax_safe=umax_safe,
                    d_umax_safe=d_umax_safe,
                    seed=args.seed,
                )
                result = cem_mpc(x0, torch.tensor(target_seq, dtype=torch.float64), li, cfg_dyn, cfg)
                u_seq = result.u_roll
                if u_seq.shape[0] < n_steps:
                    pad = u_seq[-1:].repeat(n_steps - u_seq.shape[0], 1, 1)
                    u_seq = torch.cat([u_seq, pad], dim=0)
            else:
                u_warm = torch.zeros((args.horizon, 1, 3), dtype=torch.float64)
                if currents is not None:
                    warm_raw = _slice_with_wrap(currents, s_idx, args.horizon, wrap_len=wrap_len)
                    warm_np = _clamp_currents(warm_raw, umax_safe, d_umax_safe)
                    u_warm = torch.tensor(warm_np, dtype=torch.float64).view(args.horizon, 1, 3)

                cfg = ILQRConfig(
                    horizon=args.horizon,
                    max_iter=args.max_iter,
                    w_tip=1.0,
                    w_tip_terminal=25.0,
                    w_u=5e-3,
                    w_du=5e-2,
                    max_du=d_umax_safe,
                    max_u=umax_safe,
                    linearization_backend=args.backend,
                    linearization_dense=True,
                )
                x_roll, u_roll, _, stats = run_mpc(
                    x0,
                    torch.tensor(target_seq, dtype=torch.float64),
                    li,
                    cfg_dyn,
                    cfg,
                    u_warm,
                    max_wall_sec=args.max_wall_sec,
                )
                u_seq = u_roll
                if u_seq.shape[0] < n_steps:
                    pad = u_seq[-1:].repeat(n_steps - u_seq.shape[0], 1, 1)
                    u_seq = torch.cat([u_seq, pad], dim=0)
                if any(stat.get("timed_out") for stat in stats):
                    max_wall_hit = True

            x_np, tip_np, solver_exit, residual_norm, unbounded, nonfinite, errors = _rollout_collect(
                x0,
                u_seq,
                target_seq[:-1],
                li,
                cfg_dyn,
            )
            n_collected = x_np.shape[0]
            if n_collected == 0:
                continue

            total_steps += n_collected
            fail_solver = solver_exit[:n_collected] != 0
            fail_residual = residual_norm[:n_collected] > cfg_dyn.residual_threshold
            fail_unbounded = unbounded[:n_collected]
            fail_nonfinite = nonfinite[:n_collected]
            target_slice = target_seq[:-1][:n_collected]
            fail_target = ~np.isfinite(target_slice).all(axis=1)
            fail_any = fail_solver | fail_residual | fail_unbounded | fail_nonfinite | fail_target

            if args.filter_policy == "keep_all":
                keep_mask = np.ones(n_collected, dtype=bool)
            elif args.filter_policy == "drop_unbounded":
                keep_mask = ~fail_any
            else:
                if np.any(fail_any):
                    keep_mask = np.zeros(n_collected, dtype=bool)
                    drop_reasons["policy_drop"] += int(n_collected)
                else:
                    keep_mask = np.ones(n_collected, dtype=bool)

            drop_reason_code = np.zeros(n_collected, dtype=np.int32)
            if args.filter_policy == "drop_any_fail" and np.any(fail_any):
                drop_reason_code[:] = 6
            else:
                drop_reason_code[fail_solver] = 1
                drop_reason_code[fail_residual] = 2
                drop_reason_code[fail_unbounded] = 3
                drop_reason_code[fail_nonfinite] = 4
                drop_reason_code[fail_target] = 5

            drop_mask = ~keep_mask
            drop_reasons["solver_exit"] += int((fail_solver & drop_mask).sum())
            drop_reasons["residual"] += int((fail_residual & drop_mask).sum())
            drop_reasons["unbounded"] += int((fail_unbounded & drop_mask).sum())
            drop_reasons["nonfinite"] += int((fail_nonfinite & drop_mask).sum())
            drop_reasons["target_nonfinite"] += int((fail_target & drop_mask).sum())

            all_x_raw.append(x_np[:n_collected])
            all_tip_raw.append(tip_np[:n_collected])
            all_target_raw.append(target_slice[:n_collected])
            all_u_raw.append(u_seq.detach().cpu().numpy()[:n_collected, 0])
            all_episode_raw.append(np.full((n_collected,), ep_id + episode_offset, dtype=np.int32))
            all_step_raw.append(np.arange(n_collected, dtype=np.int32))
            all_solver_raw.append(solver_exit[:n_collected])
            all_residual_raw.append(residual_norm[:n_collected])
            all_unbounded_raw.append(unbounded[:n_collected])
            all_nonfinite_raw.append(nonfinite[:n_collected])
            all_keep_raw.append(keep_mask[:n_collected])
            all_drop_reason_raw.append(drop_reason_code[:n_collected])

            kept = int(keep_mask.sum())
            dropped = int(n_collected - kept)
            kept_steps += kept
            dropped_steps += dropped

            if kept == 0:
                continue

            all_x.append(x_np[keep_mask])
            all_tip.append(tip_np[keep_mask])
            all_target.append(target_slice[keep_mask])
            all_u.append(u_seq.detach().cpu().numpy()[:n_collected, 0][keep_mask])
            all_episode.append(np.full((kept,), ep_id + episode_offset, dtype=np.int32))
            all_step.append(np.arange(n_collected, dtype=np.int32)[keep_mask])
            all_solver.append(solver_exit[:n_collected][keep_mask])
            all_residual.append(residual_norm[:n_collected][keep_mask])
            all_unbounded.append(unbounded[:n_collected][keep_mask])
            all_nonfinite.append(nonfinite[:n_collected][keep_mask])
        episode_offset += episodes

    if not all_x:
        raise SystemExit("No data collected; check indices and controller settings.")

    x_all = np.concatenate(all_x, axis=0)
    target_all = np.concatenate(all_target, axis=0)
    u_all = np.concatenate(all_u, axis=0)
    episode_all = np.concatenate(all_episode, axis=0)
    step_all = np.concatenate(all_step, axis=0)
    solver_all = np.concatenate(all_solver, axis=0)
    residual_all = np.concatenate(all_residual, axis=0)
    unbounded_all = np.concatenate(all_unbounded, axis=0)
    tip_all = np.concatenate(all_tip, axis=0)
    nonfinite_all = np.concatenate(all_nonfinite, axis=0)

    errors = np.linalg.norm(tip_all - target_all, axis=1)
    mean_error = float(np.mean(errors)) if errors.size else 0.0
    max_error = float(np.max(errors)) if errors.size else 0.0
    keep_mask_all = np.ones((x_all.shape[0],), dtype=bool)
    drop_reason_all = np.zeros((x_all.shape[0],), dtype=np.int32)

    x_raw = np.concatenate(all_x_raw, axis=0)
    target_raw = np.concatenate(all_target_raw, axis=0)
    u_raw = np.concatenate(all_u_raw, axis=0)
    episode_raw = np.concatenate(all_episode_raw, axis=0)
    step_raw = np.concatenate(all_step_raw, axis=0)
    solver_raw = np.concatenate(all_solver_raw, axis=0)
    residual_raw = np.concatenate(all_residual_raw, axis=0)
    unbounded_raw = np.concatenate(all_unbounded_raw, axis=0)
    tip_raw = np.concatenate(all_tip_raw, axis=0)
    nonfinite_raw = np.concatenate(all_nonfinite_raw, axis=0)
    keep_raw = np.concatenate(all_keep_raw, axis=0)
    drop_reason_raw = np.concatenate(all_drop_reason_raw, axis=0)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_path = args.output_npz or os.path.join("output_data", f"bc_dataset_{stamp}.npz")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    dataset_name = ",".join(dataset_names) if dataset_names else os.path.basename(args.dataset)
    np.savez(
        out_path,
        x_t=x_all,
        target_xyz=target_all,
        u_t=u_all,
        dt=dt,
        Li=li_val,
        episode_id=episode_all,
        step_id=step_all,
        solver_exit=solver_all,
        residual_norm=residual_all,
        unbounded_detected=unbounded_all,
        nonfinite_detected=nonfinite_all,
        keep_mask=keep_mask_all,
        drop_reason_code=drop_reason_all,
        x_t_all=x_raw,
        target_xyz_all=target_raw,
        u_t_all=u_raw,
        episode_id_all=episode_raw,
        step_id_all=step_raw,
        solver_exit_all=solver_raw,
        residual_norm_all=residual_raw,
        unbounded_detected_all=unbounded_raw,
        nonfinite_detected_all=nonfinite_raw,
        keep_mask_all=keep_raw,
        drop_reason_code_all=drop_reason_raw,
        dataset_name=dataset_name,
        start_idx=args.start_idx,
        end_idx=args.end_idx,
        N_total=total_n_total,
        ramp_len=ramp_len_report,
        umax=umax_safe,
        d_umax=d_umax_safe,
        n_steps_total=total_steps,
        n_steps_kept=kept_steps,
        n_steps_dropped=dropped_steps,
        drop_reasons_counts=json.dumps(drop_reasons),
        filter_policy=args.filter_policy,
    )

    report_path = args.report_path or os.path.join(
        "docs",
        "control",
        "run_reports",
        f"generate_bc_dataset_{stamp}.md",
    )
    write_run_report(
        report_path=report_path,
        command=" ".join(sys.argv),
        dataset_name=dataset_name,
        start_idx=args.start_idx,
        end_idx=args.end_idx,
        n_total=total_n_total,
        li_mm=li_report,
        dt=dt_report,
        umax=umax_safe,
        d_umax=d_umax_safe,
        mean_error=mean_error,
        max_error=max_error,
        artifact_npz=out_path,
        plot_prefix=None,
        init_from_ramp=args.init_from_ramp,
        ramp_len=ramp_len_report,
        max_wall_hit=max_wall_hit,
        extra_entries={
            "controller": args.controller,
            "filter_policy": args.filter_policy,
            "n_steps_total": total_steps,
            "n_steps_kept": kept_steps,
            "n_steps_dropped": dropped_steps,
            "drop_reasons_counts": json.dumps(drop_reasons),
        },
    )

    print(f"Saved dataset to {out_path}")
    print(
        "dataset_counts: "
        f"total={total_steps} kept={kept_steps} dropped={dropped_steps} "
        f"reasons={json.dumps(drop_reasons)}"
    )
    print(f"Run report: {report_path}")


if __name__ == "__main__":
    main()
