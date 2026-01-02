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

from crm_diffsims.control.cem_mpc import CEMConfig, cem_mpc
from crm_diffsims.control.ilqr import ILQRConfig, ilqr_solve
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


def _atomic_save_npz(path, **kwargs):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp_path = f"{path}.tmp.npz"
    np.savez(tmp_path, **kwargs)
    os.replace(tmp_path, path)


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
    return target, currents, dt, li


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


def _prepare_targets(targets, start_idx, n_steps):
    end_idx = start_idx + n_steps
    clipped = targets[start_idx:end_idx]
    if clipped.shape[0] < n_steps:
        pad = np.repeat(clipped[-1:], n_steps - clipped.shape[0], axis=0)
        clipped = np.concatenate([clipped, pad], axis=0)
    return clipped


def _rollout_prefix(x0, u_prefix, li, cfg_dyn):
    if u_prefix.shape[0] == 0:
        return x0
    x_t = x0
    for t in range(u_prefix.shape[0]):
        u_t = torch.tensor(u_prefix[t : t + 1], dtype=torch.float64).view(1, 1, 3)
        print(f"step_start mode=prefix step={t}", flush=True)
        step_start = time.time()
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        step_elapsed = time.time() - step_start
        print(f"step_end mode=prefix step={t} step_sec={step_elapsed:.2f}", flush=True)
        if step_elapsed > 2.0:
            raise RuntimeError(f"prefix step timeout at {t}")
        solver_exit = int(cache["solver_exit"].item())
        residual = float(cache["residual_norm"].item())
        if solver_exit != 0 or residual > cfg_dyn.residual_threshold:
            raise RuntimeError(f"prefix nonconvergence at {t}")
        if not torch.isfinite(x_tp1).all().item():
            raise RuntimeError(f"prefix nonfinite at {t}")
        x_t = x_tp1
    return x_t


def _rollout_open_loop_with_walltime(x0, u_seq, li, cfg_dyn, mode_name, max_wall_sec):
    x_list = []
    tip_list = []
    solver_exit = []
    residual = []
    unbounded = []
    timeout = False
    timeout_step = None
    start_time = time.time()

    x_t = x0
    for t in range(u_seq.shape[0]):
        if max_wall_sec is not None and time.time() - start_time >= max_wall_sec:
            timeout = True
            timeout_step = t
            break
        x_list.append(x_t.detach().cpu().numpy()[0])
        tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())
        u_t = torch.tensor(u_seq[t : t + 1], dtype=torch.float64).view(1, 1, 3)
        print(f"step_start mode={mode_name} step={t}", flush=True)
        step_start = time.time()
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        step_elapsed = time.time() - step_start
        print(f"step_end mode={mode_name} step={t} step_sec={step_elapsed:.2f}", flush=True)
        if step_elapsed > 2.0:
            timeout = True
            timeout_step = t
            break
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
        timeout,
        timeout_step,
    )


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
    prev_u,
    max_wall_sec,
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
    start_time = time.time()

    x_t = x0
    for t in range(targets.shape[0]):
        if max_wall_sec is not None and time.time() - start_time >= max_wall_sec:
            timeout = True
            timeout_step = t
            break
        x_list.append(x_t.detach().cpu().numpy()[0])
        tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())
        u_pred, u_used_t = _policy_step(
            model,
            x_t,
            targets[t],
            prev_u,
            umax,
            d_umax,
            include_tip,
            include_prev_u,
        )
        u_raw.append(u_pred)
        u_used.append(u_used_t)
        prev_u = u_used_t.copy()

        u_t = torch.tensor(u_used_t, dtype=torch.float64).view(1, 1, 3)
        print(f"step_start mode=policy_only step={t}", flush=True)
        step_start = time.time()
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        step_elapsed = time.time() - step_start
        print(f"step_end mode=policy_only step={t} step_sec={step_elapsed:.2f}", flush=True)
        if step_elapsed > 2.0:
            timeout = True
            timeout_step = t
            break
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


def _build_policy_warmstart(model, x_t, targets, li, cfg_dyn, umax, d_umax, include_tip, include_prev_u, u_prev_init, horizon):
    u_seq = []
    u_prev = u_prev_init.copy()
    x_pred = x_t
    for t in range(horizon):
        target_t = targets[min(t, targets.shape[0] - 1)]
        _, u_used = _policy_step(model, x_pred, target_t, u_prev, umax, d_umax, include_tip, include_prev_u)
        u_seq.append(u_used)
        u_prev = u_used.copy()
        u_t = torch.tensor(u_used, dtype=torch.float64).view(1, 1, 3)
        x_pred, cache = crm_step_v1_3_forward(x_pred, u_t, li, cfg_dyn)
        if int(cache["solver_exit"].item()) != 0:
            break
    if not u_seq:
        u_seq = [u_prev_init]
    if len(u_seq) < horizon:
        pad = np.repeat(u_seq[-1][None, :], horizon - len(u_seq), axis=0)
        u_seq = np.concatenate([np.asarray(u_seq), pad], axis=0)
    else:
        u_seq = np.asarray(u_seq)
    return torch.tensor(u_seq, dtype=torch.float64).view(horizon, 1, 3)


def _mpc_loop(x0, targets, li, cfg_dyn, cfg, u_warm_fn, mpc_steps, max_wall_sec, mode_name):
    x_roll = [x0]
    u_roll = []
    tip_list = []
    solver_exit = []
    residual = []
    unbounded = []
    timeout = False
    timeout_step = None
    start_time = time.time()

    n_steps = targets.shape[0]
    replan_interval = max(int(mpc_steps or 1), 1)
    t = 0
    while t < n_steps:
        if max_wall_sec is not None and time.time() - start_time >= max_wall_sec:
            timeout = True
            timeout_step = t
            break
        horizon_targets = targets[t : t + cfg.horizon + 1]
        if horizon_targets.shape[0] == 0:
            break
        effective_horizon = max(int(horizon_targets.shape[0] - 1), 1)
        if horizon_targets.shape[0] < cfg.horizon + 1:
            pad = horizon_targets[-1:].repeat(cfg.horizon + 1 - horizon_targets.shape[0], 1)
            horizon_targets = torch.cat([horizon_targets, pad], dim=0)
        u_init = u_warm_fn(x_roll[-1], horizon_targets)
        if u_init.shape[0] > effective_horizon:
            u_init = u_init[:effective_horizon]
        elif u_init.shape[0] < effective_horizon:
            pad = u_init[-1:].repeat(effective_horizon - u_init.shape[0], 1, 1)
            u_init = torch.cat([u_init, pad], dim=0)
        result = ilqr_solve(
            x_roll[-1],
            u_init,
            horizon_targets,
            li,
            cfg_dyn,
            cfg,
            max_wall_sec=max_wall_sec,
            start_time=time.time() if max_wall_sec is not None else None,
        )
        if result.timed_out:
            timeout = True
            timeout_step = t
            break
        u_seq = result.u_seq
        n_apply = min(replan_interval, u_seq.shape[0], n_steps - t)
        if n_apply == 0:
            break
        for k in range(n_apply):
            if max_wall_sec is not None and time.time() - start_time >= max_wall_sec:
                timeout = True
                timeout_step = t
                break
            u_apply = u_seq[k : k + 1]
            u_roll.append(u_apply)
            tip_list.append(x_roll[-1][0, TIP_IDX].detach().cpu().numpy())
            print(f"step_start mode={mode_name} step={t}", flush=True)
            step_start = time.time()
            x_next, cache = crm_step_v1_3_forward(x_roll[-1], u_apply, li, cfg_dyn)
            step_elapsed = time.time() - step_start
            print(f"step_end mode={mode_name} step={t} step_sec={step_elapsed:.2f}", flush=True)
            if step_elapsed > 2.0:
                timeout = True
                timeout_step = t
                break
            exit_flag = int(cache["solver_exit"].item())
            res = float(cache["residual_norm"].item())
            solver_exit.append(exit_flag)
            residual.append(res)
            unbounded_flag = not torch.isfinite(x_next).all().item()
            unbounded.append(unbounded_flag)
            if exit_flag != 0 or res > cfg_dyn.residual_threshold:
                t = n_steps
                break
            if unbounded_flag:
                t = n_steps
                break
            x_roll.append(x_next)
            t += 1
        if timeout:
            break

    return (
        torch.cat(x_roll, dim=0).detach().cpu().numpy(),
        np.asarray(tip_list, dtype=np.float64),
        torch.cat(u_roll, dim=0).detach().cpu().numpy() if u_roll else np.zeros((0, 1, 3)),
        np.asarray(solver_exit, dtype=np.int32),
        np.asarray(residual, dtype=np.float64),
        np.asarray(unbounded, dtype=bool),
        timeout,
        timeout_step,
    )


def _keep_mask(tip, target, solver_exit, residual, unbounded, residual_threshold):
    lengths = [
        tip.shape[0],
        target.shape[0],
        solver_exit.shape[0],
        residual.shape[0],
        unbounded.shape[0],
    ]
    n = min(lengths) if lengths else 0
    if n == 0:
        return np.zeros((0,), dtype=bool)
    tip = tip[:n]
    target = target[:n]
    solver_exit = solver_exit[:n]
    residual = residual[:n]
    unbounded = unbounded[:n]
    nonfinite = ~np.isfinite(tip).all(axis=1)
    target_nonfinite = ~np.isfinite(target).all(axis=1)
    fail = (
        unbounded.astype(bool)
        | nonfinite
        | target_nonfinite
        | (solver_exit != 0)
        | (residual > residual_threshold)
    )
    return ~fail


def _error_metrics(tip, target, keep):
    if keep.sum() == 0:
        return float("nan"), float("nan"), np.zeros((0,))
    err = np.linalg.norm(tip[keep] - target[: tip.shape[0]][keep], axis=1)
    return float(np.mean(err)), float(np.max(err)), err


def _solver_exit_counts(solver_exit):
    if solver_exit.size == 0:
        return {}
    unique, counts = np.unique(solver_exit, return_counts=True)
    return {int(k): int(v) for k, v in zip(unique, counts)}


def _normalize_u(u_seq):
    if u_seq is None:
        return np.zeros((0, 3), dtype=np.float64)
    if u_seq.ndim == 3 and u_seq.shape[1] == 1:
        return u_seq[:, 0]
    return u_seq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, default="data/dyn_fk_ramp_circle1_hold1.npz")
    parser.add_argument("--start-idx", type=int, default=0)
    parser.add_argument("--window-len", type=int, default=50)
    parser.add_argument("--model", type=str, default=None)
    parser.add_argument(
        "--modes",
        type=str,
        default="replay_only,policy_only,policy_plus_mpc,cem_mpc",
    )
    parser.add_argument("--preroll-steps", type=int, default=0)
    parser.add_argument("--max-wall-sec", type=float, default=60.0)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    if args.window_len <= 0:
        raise SystemExit("--window-len must be > 0")

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    ext_dir = os.path.join(repo_root, "build_torch_ext_dyn")
    os.environ["TORCH_EXTENSIONS_DIR"] = ext_dir
    print(f"TORCH_EXTENSIONS_DIR={ext_dir}", flush=True)
    lock_path = os.path.join(ext_dir, "lock")
    if os.path.exists(lock_path):
        print(f"Lock exists: {lock_path}", flush=True)
        raise SystemExit(2)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    targets, currents, dt, li_val = _load_dataset(args.dataset)
    if currents is None:
        raise SystemExit("dataset currents required")
    umax_safe, d_umax_safe = _load_safe_bounds()

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0_default, _ = build_state()
    li = torch.tensor([li_val], dtype=torch.float64)

    base_start_idx = args.start_idx
    preroll_steps = max(int(args.preroll_steps), 0)
    eval_start_idx = base_start_idx + preroll_steps
    eval_end_idx = eval_start_idx + int(args.window_len)

    x0 = x0_default.clone()
    if base_start_idx > 0:
        prefix_u = _clamp_currents(currents[:base_start_idx], umax_safe, d_umax_safe)
        x0 = _rollout_prefix(x0, prefix_u, li, cfg_dyn)
    if preroll_steps > 0:
        preroll_u = _clamp_currents(
            currents[base_start_idx : base_start_idx + preroll_steps],
            umax_safe,
            d_umax_safe,
        )
        x0 = _rollout_prefix(x0, preroll_u, li, cfg_dyn)
    else:
        preroll_u = np.zeros((0, 3), dtype=np.float64)

    target_seq = _prepare_targets(targets, eval_start_idx, args.window_len)

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]

    model = None
    model_meta = {}
    if any(m in modes for m in ("policy_only", "policy_plus_mpc")):
        if args.model is None:
            if os.path.exists(os.path.join("output_data", "models")):
                candidates = [
                    os.path.join("output_data", "models", name)
                    for name in os.listdir(os.path.join("output_data", "models"))
                    if name.endswith(".pt")
                ]
                if candidates:
                    args.model = max(candidates, key=os.path.getmtime)
        if args.model and os.path.exists(args.model):
            model, model_meta = _load_model(args.model)
        else:
            print("No model available; skipping policy modes.", flush=True)
            modes = [m for m in modes if m not in ("policy_only", "policy_plus_mpc")]

    include_tip = bool(model_meta.get("include_tip", False)) if model_meta else False
    include_prev_u = bool(model_meta.get("include_prev_u", False)) if model_meta else False
    prev_u_init = preroll_u[-1].copy() if preroll_u.shape[0] else np.zeros((3,), dtype=np.float64)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    for mode in modes:
        start_time = time.time()
        if mode == "replay_only":
            u_raw = currents[eval_start_idx:eval_end_idx]
            u_seq = _clamp_currents(u_raw, umax_safe, d_umax_safe)
            x_np, tip_np, solver_exit, residual, unbounded, timeout, timeout_step = _rollout_open_loop_with_walltime(
                x0, u_seq, li, cfg_dyn, "replay_only", args.max_wall_sec
            )
            u_used = u_seq[: tip_np.shape[0]]
            u_raw_used = u_used.copy()
        elif mode == "policy_only":
            if model is None:
                continue
            (
                x_np,
                tip_np,
                u_used,
                u_raw_used,
                solver_exit,
                residual,
                unbounded,
                timeout,
                timeout_step,
            ) = _rollout_policy_only(
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
                args.max_wall_sec,
            )
        elif mode == "policy_plus_mpc":
            if model is None:
                continue
            cfg = ILQRConfig(
                horizon=min(10, int(args.window_len)),
                max_iter=1,
                w_tip=1.0,
                w_tip_terminal=10.0,
                w_u=1e-2,
                w_du=1e-1,
                max_du=d_umax_safe,
                max_u=umax_safe,
                linearization_backend="v1_3",
                linearization_dense=True,
            )

            def _policy_warm(x_t, horizon_targets):
                return _build_policy_warmstart(
                    model,
                    x_t,
                    horizon_targets[:, :],
                    li,
                    cfg_dyn,
                    umax_safe,
                    d_umax_safe,
                    include_tip,
                    include_prev_u,
                    prev_u_init,
                    cfg.horizon,
                )

            (
                x_np,
                tip_np,
                u_used,
                solver_exit,
                residual,
                unbounded,
                timeout,
                timeout_step,
            ) = _mpc_loop(
                x0,
                torch.tensor(target_seq, dtype=torch.float64),
                li,
                cfg_dyn,
                cfg,
                _policy_warm,
                mpc_steps=1,
                max_wall_sec=args.max_wall_sec,
                mode_name="policy_plus_mpc",
            )
            u_raw_used = u_used.copy()
        elif mode == "cem_mpc":
            cfg = CEMConfig(
                horizon=min(10, int(args.window_len)),
                mpc_steps=int(args.window_len),
                num_samples=32,
                num_elites=6,
                cem_iters=2,
                u_std_init=0.02,
                umax_safe=umax_safe,
                d_umax_safe=d_umax_safe,
                seed=args.seed,
            )
            result = cem_mpc(x0, torch.tensor(target_seq, dtype=torch.float64), li, cfg_dyn, cfg)
            u_seq = result.u_roll.detach().cpu().numpy()[:, 0]
            if u_seq.shape[0] < args.window_len:
                pad = np.repeat(u_seq[-1:], args.window_len - u_seq.shape[0], axis=0)
                u_seq = np.concatenate([u_seq, pad], axis=0)
            u_seq = _clamp_currents(u_seq, umax_safe, d_umax_safe)
            x_np, tip_np, solver_exit, residual, unbounded, timeout, timeout_step = _rollout_open_loop_with_walltime(
                x0, u_seq, li, cfg_dyn, "cem_mpc", args.max_wall_sec
            )
            u_used = u_seq[: tip_np.shape[0]]
            u_raw_used = u_used.copy()
        else:
            print(f"Unknown mode: {mode}", flush=True)
            continue

        u_used = _normalize_u(u_used)
        u_raw_used = _normalize_u(u_raw_used)

        keep = _keep_mask(tip_np, target_seq, solver_exit, residual, unbounded, cfg_dyn.residual_threshold)
        mean_err, max_err, errors = _error_metrics(tip_np, target_seq, keep)
        runtime_sec = time.time() - start_time
        solver_counts = _solver_exit_counts(solver_exit)
        unbounded_count = int(unbounded.sum()) if unbounded.size else 0
        kept_steps = int(keep.sum())
        attempted = int(tip_np.shape[0])

        stamp = f"{ts}_{mode}"
        out_npz = os.path.join("output_data", f"scoreboard_{stamp}.npz")
        report_path = os.path.join("docs", "control", "run_reports", f"scoreboard_{stamp}.md")
        _atomic_save_npz(
            out_npz,
            x_t=x_np,
            tip_xyz=tip_np,
            target_xyz=target_seq[: tip_np.shape[0]],
            u_t=u_used,
            u_raw=u_raw_used,
            solver_exit=solver_exit,
            residual_norm=residual,
            unbounded_detected=unbounded,
            keep_mask=keep,
            error_mm=errors,
            dt=float(dt),
            Li=float(li_val),
            umax=float(umax_safe),
            d_umax=float(d_umax_safe),
        )

        write_run_report(
            report_path=report_path,
            command=" ".join(sys.argv),
            dataset_name=os.path.basename(args.dataset),
            start_idx=eval_start_idx,
            end_idx=eval_end_idx,
            n_total=int(targets.shape[0]),
            li_mm=float(li_val),
            dt=float(dt),
            umax=float(umax_safe),
            d_umax=float(d_umax_safe),
            mean_error=mean_err,
            max_error=max_err,
            artifact_npz=out_npz,
            plot_prefix=None,
            init_from_ramp=False,
            ramp_len=0,
            max_wall_hit=bool(timeout),
            extra_entries={
                "mode": mode,
                "kept_steps": kept_steps,
                "attempted_steps": attempted,
                "unbounded_count": unbounded_count,
                "solver_exit_counts": json.dumps(solver_counts),
                "runtime_sec": f"{runtime_sec:.2f}",
                "timeout_step": timeout_step if timeout else None,
            },
        )

        print(
            f"summary mode={mode} kept={kept_steps}/{attempted} "
            f"mean_mm={mean_err:.3f} max_mm={max_err:.3f} "
            f"unbounded={unbounded_count} solver_exit_counts={json.dumps(solver_counts)} "
            f"runtime_sec={runtime_sec:.2f}",
            flush=True,
        )
        print(f"artifact_npz={out_npz}", flush=True)
        print(f"report_path={report_path}", flush=True)


if __name__ == "__main__":
    main()
