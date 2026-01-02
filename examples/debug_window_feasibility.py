import argparse
import json
import os
import sys
import time
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

from crm_diffsims.control.run_report import write_run_report
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


TIP_IDX = slice(24, 27)


class MLPPolicy(torch.nn.Module):
    def __init__(self, input_dim, hidden_sizes):
        super().__init__()
        layers = []
        dim = input_dim
        for h in hidden_sizes:
            layers.append(torch.nn.Linear(dim, h))
            layers.append(torch.nn.ReLU())
            dim = h
        layers.append(torch.nn.Linear(dim, 3))
        self.net = torch.nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def _load_model(path):
    payload = torch.load(path, map_location="cpu")
    hidden_sizes = payload["hidden_sizes"]
    input_dim = payload["input_dim"]
    model = MLPPolicy(input_dim, hidden_sizes)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload


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


def _load_safe_bounds():
    path = os.path.join("docs", "control", "safe_bounds.json")
    if not os.path.exists(path):
        return 0.1, 0.01
    with open(path, "r", encoding="ascii") as f:
        payload = json.load(f)
    return float(payload.get("umax_safe", 0.1)), float(payload.get("d_umax_safe", 0.01))


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


def _rollout_open_loop(x0, u_seq, li, cfg_dyn, max_steps=None, max_wall_sec=None, start_time=None):
    x_list = []
    tip_list = []
    solver_exit = []
    residual = []
    unbounded = []
    timeout = False
    timeout_step = None

    x_t = x0
    n_steps = u_seq.shape[0]
    if max_steps is not None:
        n_steps = min(n_steps, int(max_steps))
    for t in range(n_steps):
        if max_wall_sec is not None and start_time is not None:
            if time.time() - start_time >= max_wall_sec:
                timeout = True
                timeout_step = t
                break
        x_list.append(x_t.detach().cpu().numpy()[0])
        tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())
        u_t = torch.tensor(u_seq[t : t + 1], dtype=torch.float64).view(1, 1, 3)
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        exit_flag = int(cache["solver_exit"].item())
        res = float(cache["residual_norm"].item())
        solver_exit.append(exit_flag)
        residual.append(res)
        unbounded_flag = not torch.isfinite(x_tp1).all().item()
        unbounded.append(unbounded_flag)

        if exit_flag != 0 or res > cfg_dyn.residual_threshold:
            break
        if unbounded_flag:
            break
        x_t = x_tp1

    return (
        np.asarray(x_list, dtype=np.float64),
        np.asarray(tip_list, dtype=np.float64),
        np.asarray(solver_exit, dtype=np.int32),
        np.asarray(residual, dtype=np.float64),
        np.asarray(unbounded, dtype=bool),
        x_t,
        timeout,
        timeout_step,
    )


def _prepare_targets(targets, start_idx, n_steps):
    end_idx = start_idx + n_steps
    clipped = targets[start_idx:end_idx]
    if clipped.shape[0] < n_steps:
        pad = np.repeat(clipped[-1:], n_steps - clipped.shape[0], axis=0)
        clipped = np.concatenate([clipped, pad], axis=0)
    return clipped


def _policy_step(model, x_t, target_t, u_prev, umax, d_umax, include_tip, include_prev_u):
    if include_tip:
        state = x_t[0, TIP_IDX].detach().cpu().float()
    else:
        state = x_t[0].detach().cpu().float()
    inputs = [state, torch.tensor(target_t, dtype=torch.float32)]
    if include_prev_u:
        inputs.append(torch.tensor(u_prev, dtype=torch.float32))
    x_in = torch.cat(inputs, dim=0).unsqueeze(0)
    with torch.no_grad():
        u_pred = model(x_in)[0].cpu().numpy()
    u_clamped = np.clip(u_pred, -umax, umax)
    delta = u_clamped - u_prev
    delta = np.clip(delta, -d_umax, d_umax)
    u_used = u_prev + delta
    return u_pred, u_used


def _rollout_policy_only(
    model,
    x0,
    targets,
    li,
    cfg_dyn,
    umax,
    d_umax,
    include_tip,
    include_prev_u,
    prev_u_init,
    max_steps=None,
    max_wall_sec=None,
    start_time=None,
):
    x_list = []
    tip_list = []
    u_used = []
    u_raw = []
    solver_exit = []
    residual = []
    unbounded = []
    timeout = False
    timeout_step = None

    x_t = x0
    u_prev = prev_u_init.copy()
    n_steps = targets.shape[0]
    if max_steps is not None:
        n_steps = min(n_steps, int(max_steps))
    for t in range(n_steps):
        if max_wall_sec is not None and start_time is not None:
            if time.time() - start_time >= max_wall_sec:
                timeout = True
                timeout_step = t
                break
        x_list.append(x_t.detach().cpu().numpy()[0])
        tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())
        u_pred, u_used_t = _policy_step(
            model,
            x_t,
            targets[t],
            u_prev,
            umax,
            d_umax,
            include_tip,
            include_prev_u,
        )
        u_raw.append(u_pred)
        u_used.append(u_used_t)
        u_prev = u_used_t.copy()

        u_t = torch.tensor(u_used_t, dtype=torch.float64).view(1, 1, 3)
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        exit_flag = int(cache["solver_exit"].item())
        res = float(cache["residual_norm"].item())
        solver_exit.append(exit_flag)
        residual.append(res)
        unbounded_flag = not torch.isfinite(x_tp1).all().item()
        unbounded.append(unbounded_flag)

        if exit_flag != 0 or res > cfg_dyn.residual_threshold:
            break
        if unbounded_flag:
            break
        x_t = x_tp1

    return (
        np.asarray(x_list, dtype=np.float64),
        np.asarray(tip_list, dtype=np.float64),
        np.asarray(u_used, dtype=np.float64),
        np.asarray(u_raw, dtype=np.float64),
        np.asarray(solver_exit, dtype=np.int32),
        np.asarray(residual, dtype=np.float64),
        np.asarray(unbounded, dtype=bool),
        timeout,
        timeout_step,
    )


def _first_failure_reason(tip, target, solver_exit, residual, unbounded, residual_threshold):
    n_steps = tip.shape[0]
    for idx in range(n_steps):
        if unbounded[idx]:
            return idx, "unbounded"
        if not np.isfinite(tip[idx]).all():
            return idx, "nonfinite_tip"
        if not np.isfinite(target[idx]).all():
            return idx, "target_nonfinite"
        if solver_exit[idx] != 0:
            return idx, "solver_exit"
        if residual[idx] > residual_threshold:
            return idx, "residual"
    return None, "none"


def _summarize_run(label, tip, target, solver_exit, residual, unbounded, residual_threshold):
    total = tip.shape[0]
    if total == 0:
        stats = {
            "total_steps": 0,
            "kept_steps": 0,
            "solver_exit": 0,
            "residual": 0,
            "unbounded": 0,
            "nonfinite_tip": 0,
            "target_nonfinite": 0,
            "first_fail_step": None,
            "first_fail_reason": "no_steps",
        }
        print(f"{label}: no steps")
        return stats

    nonfinite_tip = ~np.isfinite(tip).all(axis=1)
    target_nonfinite = ~np.isfinite(target[:total]).all(axis=1)
    fail = (
        unbounded.astype(bool)
        | nonfinite_tip
        | target_nonfinite
        | (solver_exit != 0)
        | (residual > residual_threshold)
    )
    kept = int((~fail).sum())
    first_step, first_reason = _first_failure_reason(
        tip,
        target[:total],
        solver_exit,
        residual,
        unbounded.astype(bool),
        residual_threshold,
    )

    stats = {
        "total_steps": int(total),
        "kept_steps": kept,
        "solver_exit": int((solver_exit != 0).sum()),
        "residual": int((residual > residual_threshold).sum()),
        "unbounded": int(unbounded.astype(bool).sum()),
        "nonfinite_tip": int(nonfinite_tip.sum()),
        "target_nonfinite": int(target_nonfinite.sum()),
        "first_fail_step": first_step,
        "first_fail_reason": first_reason,
    }
    print(f"{label}: total={stats['total_steps']} kept={stats['kept_steps']}")
    print(
        f"{label}: solver_exit={stats['solver_exit']} residual={stats['residual']} "
        f"unbounded={stats['unbounded']} nonfinite_tip={stats['nonfinite_tip']} "
        f"target_nonfinite={stats['target_nonfinite']}"
    )
    print(f"{label}: first_fail_step={stats['first_fail_step']} reason={stats['first_fail_reason']}")
    return stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, required=True)
    parser.add_argument("--start-idx", type=int, required=True)
    parser.add_argument("--end-idx", type=int, required=True)
    parser.add_argument("--preroll-steps", type=int, default=0)
    parser.add_argument("--max-steps", type=int, default=80)
    parser.add_argument("--max-wall-sec", type=float, default=60.0)
    parser.add_argument("--init-mode", choices=["dataset_x0", "init_from_ramp", "zero"], default="dataset_x0")
    parser.add_argument("--ramp-len", type=int, default=None)
    parser.add_argument("--model", type=str, default=None)
    parser.add_argument("--output-prefix", type=str, default=None)
    parser.add_argument("--report-path", type=str, default=None)
    args = parser.parse_args()

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0_default, _ = build_state()

    if args.model is None:
        model_dir = os.path.join("output_data", "models")
        if not os.path.isdir(model_dir):
            raise SystemExit("--model is required (no models found)")
        candidates = [os.path.join(model_dir, name) for name in os.listdir(model_dir) if name.endswith(".pt")]
        if not candidates:
            raise SystemExit("--model is required (no .pt models found)")
        args.model = max(candidates, key=os.path.getmtime)
        print(f"Using latest model: {args.model}")

    model, model_meta = _load_model(args.model)
    include_tip = bool(model_meta.get("include_tip", False))
    include_prev_u = bool(model_meta.get("include_prev_u", False))

    targets, currents, dt, li_val, segments = _load_dataset(args.dataset)
    if currents is None:
        raise SystemExit("dataset currents required")

    li = torch.tensor([li_val], dtype=torch.float64)
    umax_safe, d_umax_safe = _load_safe_bounds()

    ramp_len = _infer_ramp_len(segments, args.ramp_len) if args.init_mode == "init_from_ramp" else 0
    base_start_idx = args.start_idx
    if args.init_mode == "init_from_ramp" and ramp_len > 0:
        base_start_idx = max(base_start_idx, ramp_len)

    x0 = x0_default.clone()
    start_time = time.time()
    if args.init_mode == "init_from_ramp" and ramp_len > 0:
        ramp_u = _clamp_currents(currents[:ramp_len], umax_safe, d_umax_safe)
        _, _, _, _, _, x0, _, _ = _rollout_open_loop(
            x0,
            ramp_u,
            li,
            cfg_dyn,
            max_steps=args.max_steps,
            max_wall_sec=args.max_wall_sec,
            start_time=start_time,
        )
        if base_start_idx > ramp_len:
            prefix_u = _clamp_currents(currents[ramp_len:base_start_idx], umax_safe, d_umax_safe)
            _, _, _, _, _, x0, _, _ = _rollout_open_loop(
                x0,
                prefix_u,
                li,
                cfg_dyn,
                max_steps=args.max_steps,
                max_wall_sec=args.max_wall_sec,
                start_time=start_time,
            )
    elif args.init_mode == "dataset_x0" and base_start_idx > 0:
        prefix_u = _clamp_currents(currents[:base_start_idx], umax_safe, d_umax_safe)
        _, _, _, _, _, x0, _, _ = _rollout_open_loop(
            x0,
            prefix_u,
            li,
            cfg_dyn,
            max_steps=args.max_steps,
            max_wall_sec=args.max_wall_sec,
            start_time=start_time,
        )

    preroll_steps = max(int(args.preroll_steps), 0)
    preroll_start_idx = base_start_idx
    preroll_end_idx = base_start_idx + preroll_steps
    preroll_u = _clamp_currents(currents[preroll_start_idx:preroll_end_idx], umax_safe, d_umax_safe)
    preroll_x, preroll_tip, preroll_solver, preroll_residual, preroll_unbounded, x0, preroll_timeout, preroll_timeout_step = _rollout_open_loop(
        x0,
        preroll_u,
        li,
        cfg_dyn,
        max_steps=args.max_steps,
        max_wall_sec=args.max_wall_sec,
        start_time=start_time,
    )

    window = max(args.end_idx - base_start_idx, 0)
    start_idx = base_start_idx + preroll_steps
    end_idx = start_idx + window
    target_seq = _prepare_targets(targets, start_idx, window)
    expert_u = _clamp_currents(currents[start_idx : start_idx + window], umax_safe, d_umax_safe)
    if preroll_u.shape[0]:
        prev_u_init = preroll_u[-1].copy()
    else:
        prev_u_init = expert_u[0] if expert_u.shape[0] else np.zeros((3,), dtype=np.float64)

    exp = _rollout_open_loop(
        x0,
        expert_u,
        li,
        cfg_dyn,
        max_steps=args.max_steps,
        max_wall_sec=args.max_wall_sec,
        start_time=start_time,
    )
    pol = _rollout_policy_only(
        model,
        x0,
        target_seq,
        li,
        cfg_dyn,
        umax_safe,
        d_umax_safe,
        include_tip,
        include_prev_u,
        prev_u_init,
        max_steps=args.max_steps,
        max_wall_sec=args.max_wall_sec,
        start_time=start_time,
    )

    expert_stats = _summarize_run(
        "expert",
        exp[1],
        target_seq,
        exp[2],
        exp[3],
        exp[4],
        cfg_dyn.residual_threshold,
    )
    policy_stats = _summarize_run(
        "policy",
        pol[1],
        target_seq,
        pol[4],
        pol[5],
        pol[6],
        cfg_dyn.residual_threshold,
    )

    err_policy = (
        np.linalg.norm(pol[1] - target_seq[: pol[1].shape[0]], axis=1) if pol[1].shape[0] else np.zeros((0,))
    )

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_prefix = args.output_prefix or f"debug_window_feasibility_{stamp}"
    out_path = os.path.join("output_data", f"{out_prefix}.npz")

    np.savez(
        out_path,
        dataset_name=os.path.basename(args.dataset),
        start_idx=start_idx,
        end_idx=end_idx,
        base_start_idx=base_start_idx,
        preroll_steps=preroll_steps,
        preroll_start_idx=preroll_start_idx,
        preroll_end_idx=preroll_end_idx,
        target_xyz=target_seq,
        x_preroll=preroll_x,
        tip_preroll=preroll_tip,
        u_preroll=preroll_u,
        solver_exit_preroll=preroll_solver,
        residual_preroll=preroll_residual,
        unbounded_preroll=preroll_unbounded,
        timeout_preroll=preroll_timeout,
        timeout_step_preroll=preroll_timeout_step,
        x_expert=exp[0],
        tip_expert=exp[1],
        u_expert=expert_u[: exp[1].shape[0]],
        solver_exit_expert=exp[2],
        residual_expert=exp[3],
        unbounded_expert=exp[4],
        timeout_expert=exp[6],
        timeout_step_expert=exp[7],
        x_policy=pol[0],
        tip_policy=pol[1],
        u_policy=pol[2],
        u_raw_policy=pol[3],
        solver_exit_policy=pol[4],
        residual_policy=pol[5],
        unbounded_policy=pol[6],
        timeout_policy=pol[7],
        timeout_step_policy=pol[8],
        error_policy=err_policy,
        dt=dt,
        Li=li_val,
        init_mode=args.init_mode,
        model_path=args.model,
        umax=umax_safe,
        d_umax=d_umax_safe,
        max_steps=args.max_steps,
        max_wall_sec=args.max_wall_sec,
    )

    report_path = args.report_path or os.path.join(
        "docs",
        "control",
        "run_reports",
        f"debug_window_feasibility_{stamp}.md",
    )

    write_run_report(
        report_path=report_path,
        command=" ".join(sys.argv),
        dataset_name=os.path.basename(args.dataset),
        start_idx=start_idx,
        end_idx=end_idx,
        n_total=targets.shape[0],
        li_mm=li_val,
        dt=dt,
        umax=umax_safe,
        d_umax=d_umax_safe,
        mean_error=float(np.mean(err_policy)) if err_policy.size else float("nan"),
        max_error=float(np.max(err_policy)) if err_policy.size else float("nan"),
        artifact_npz=out_path,
        plot_prefix=None,
        init_from_ramp=args.init_mode == "init_from_ramp",
        ramp_len=ramp_len,
        max_wall_hit=False,
        extra_entries={
            "model_path": args.model,
            "base_start_idx": base_start_idx,
            "preroll_steps": preroll_steps,
            "preroll_start_idx": preroll_start_idx,
            "preroll_end_idx": preroll_end_idx,
            "max_steps": args.max_steps,
            "max_wall_sec": args.max_wall_sec,
            "timeout_preroll": preroll_timeout,
            "timeout_preroll_step": preroll_timeout_step,
            "timeout_expert": exp[6],
            "timeout_expert_step": exp[7],
            "timeout_policy": pol[7],
            "timeout_policy_step": pol[8],
            "expert_stats": expert_stats,
            "policy_stats": policy_stats,
        },
    )

    print(f"Saved {out_path}")
    print(f"Report {report_path}")


if __name__ == "__main__":
    main()
