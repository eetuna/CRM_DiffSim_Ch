import argparse
import json
import os
import signal
import sys
import time
import traceback
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

from crm_diffsims.control.ilqr import ILQRConfig, ilqr_solve
from crm_diffsims.control.run_report import write_run_report
from crm_diffsims.control import plot_utils
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


def _atomic_write_report(*, report_path, **kwargs):
    tmp_path = f"{report_path}.tmp"
    write_run_report(report_path=tmp_path, **kwargs)
    os.replace(tmp_path, report_path)


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


def _load_safe_bounds_optional():
    path = os.path.join("docs", "control", "safe_bounds.json")
    if not os.path.exists(path):
        return None, None
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


def _rollout_prefix(x0, u_prefix, li, cfg_dyn, mode_name="prefix"):
    if u_prefix.shape[0] == 0:
        return x0, False, None
    x_t = x0
    for t in range(u_prefix.shape[0]):
        u_t = torch.tensor(u_prefix[t : t + 1], dtype=torch.float64).view(1, 1, 3)
        print(f"step_start mode={mode_name} step={t}", flush=True)
        step_start = time.time()
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        step_elapsed = time.time() - step_start
        print(f"step_end mode={mode_name} step={t} step_sec={step_elapsed:.2f}", flush=True)
        if time.time() - step_start > 2.0:
            return x_t, True, t
        solver_exit = int(cache["solver_exit"].item())
        residual = float(cache["residual_norm"].item())
        if solver_exit != 0 or residual > cfg_dyn.residual_threshold:
            raise RuntimeError(f"nonconvergence during prefix at step {t}")
        if not torch.isfinite(x_tp1).all().item():
            raise RuntimeError(f"nonfinite during prefix at step {t}")
        x_t = x_tp1
    return x_t, False, None


def _rollout_open_loop(x0, u_seq, li, cfg_dyn, mode_name="open_loop"):
    x_list = []
    tip_list = []
    solver_exit = []
    residual = []
    unbounded = []

    x_t = x0
    timeout = False
    timeout_step = None
    for t in range(u_seq.shape[0]):
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


def _rollout_policy_only(model, x0, targets, li, cfg_dyn, umax, d_umax, include_tip, include_prev_u, prev_u_init):
    x_list = []
    tip_list = []
    u_used = []
    u_raw = []
    solver_exit = []
    residual = []
    unbounded = []

    x_t = x0
    u_prev = prev_u_init.copy()
    return _rollout_policy_only_with_walltime(
        model,
        x_t,
        targets,
        li,
        cfg_dyn,
        umax,
        d_umax,
        include_tip,
        include_prev_u,
        u_prev,
        None,
        None,
    )


def _rollout_policy_only_with_walltime(
    model,
    x_t,
    targets,
    li,
    cfg_dyn,
    umax,
    d_umax,
    include_tip,
    include_prev_u,
    u_prev,
    max_wall_sec,
    start_time,
    heartbeat_steps,
    mode_name,
    rollout_tracker,
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
    timeout_reason = "none"

    kept_steps = 0
    for t in range(targets.shape[0]):
        if max_wall_sec is not None and start_time is not None:
            if time.time() - start_time >= max_wall_sec:
                timeout = True
                timeout_step = t
                timeout_reason = "max_wall_sec"
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
        print(
            f"step_start mode={mode_name} step={t} elapsed_sec={time.time() - start_time:.2f}",
            flush=True,
        )
        step_start = time.time()
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        step_elapsed = time.time() - step_start
        if rollout_tracker is not None:
            rollout_tracker[0] = time.time() - start_time
        print(
            f"step_end mode={mode_name} step={t} step_sec={step_elapsed:.2f} elapsed_sec={time.time() - start_time:.2f}",
            flush=True,
        )
        if step_elapsed > 2.0:
            timeout = True
            timeout_step = t
            timeout_reason = "step_timeout"
            break
        exit_flag = int(cache["solver_exit"].item())
        res = float(cache["residual_norm"].item())
        solver_exit.append(exit_flag)
        residual.append(res)
        unbounded_flag = not torch.isfinite(x_tp1).all().item()
        unbounded.append(unbounded_flag)

        if heartbeat_steps and (t % heartbeat_steps == 0):
            print(
                f"heartbeat mode={mode_name} step={t} elapsed_sec={time.time() - start_time:.2f} "
                f"solver_exit={exit_flag} residual={res:.3e} kept_steps={kept_steps}",
                flush=True,
            )

        if exit_flag != 0 or res > cfg_dyn.residual_threshold:
            break
        if unbounded_flag:
            break
        kept_steps += 1
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
        timeout_reason,
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
        print(f"step_start mode=policy_warmstart step={t}", flush=True)
        step_start = time.time()
        x_pred, cache = crm_step_v1_3_forward(x_pred, u_t, li, cfg_dyn)
        step_elapsed = time.time() - step_start
        print(f"step_end mode=policy_warmstart step={t} step_sec={step_elapsed:.2f}", flush=True)
        if step_elapsed > 2.0:
            break
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


def _mpc_loop(
    x0,
    targets,
    li,
    cfg_dyn,
    cfg,
    u_warm_fn,
    mpc_steps,
    max_wall_sec,
    start_time,
    heartbeat_steps,
    mode_name,
    rollout_tracker,
):
    x_roll = [x0]
    u_roll = []
    tip_list = []
    solver_exit = []
    residual = []
    unbounded = []
    timeout = False
    timeout_step = None
    timeout_reason = "none"
    exception_type = ""
    exception_msg = ""
    exception_tb = ""
    keep_flags = []
    exception_hit = False

    n_steps = targets.shape[0]
    replan_interval = max(int(mpc_steps or 1), 1)
    t = 0
    kept_steps = 0
    while t < n_steps:
        if max_wall_sec is not None and start_time is not None:
            elapsed = time.time() - start_time
            if elapsed >= max_wall_sec:
                timeout = True
                timeout_step = t
                timeout_reason = "max_wall_sec"
                break
            remaining = max(max_wall_sec - elapsed, 0.0)
        else:
            remaining = None

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
        try:
            result = ilqr_solve(
                x_roll[-1],
                u_init,
                horizon_targets,
                li,
                cfg_dyn,
                cfg,
                max_wall_sec=remaining,
                start_time=time.time() if remaining is not None else None,
            )
        except Exception as exc:
            print(f"iLQR exception in {mode_name} at step {t}: {exc}", flush=True)
            exception_type = type(exc).__name__
            exception_msg = str(exc)
            exception_tb = "\n".join(traceback.format_exc().splitlines()[:20])
            break
        if result.timed_out:
            timeout = True
            timeout_step = t
            timeout_reason = "max_wall_sec"
            break
        u_seq = result.u_seq
        n_apply = min(replan_interval, u_seq.shape[0], n_steps - t)
        if n_apply == 0:
            break
        for k in range(n_apply):
            if max_wall_sec is not None and start_time is not None:
                if time.time() - start_time >= max_wall_sec:
                    timeout = True
                    timeout_step = t
                    timeout_reason = "max_wall_sec"
                    break
            u_apply = u_seq[k : k + 1]
            u_roll.append(u_apply)
            tip_list.append(x_roll[-1][0, TIP_IDX].detach().cpu().numpy())

            print(
                f"step_start mode={mode_name} step={t} elapsed_sec={time.time() - start_time:.2f}",
                flush=True,
            )
            step_start = time.time()
            x_next, cache = crm_step_v1_3_forward(x_roll[-1], u_apply, li, cfg_dyn)
            step_elapsed = time.time() - step_start
            if rollout_tracker is not None:
                rollout_tracker[0] = time.time() - start_time
            print(
                f"step_end mode={mode_name} step={t} step_sec={step_elapsed:.2f} elapsed_sec={time.time() - start_time:.2f}",
                flush=True,
            )
            if step_elapsed > 2.0:
                timeout = True
                timeout_step = t
                timeout_reason = "step_timeout"
                break

            exit_flag = int(cache["solver_exit"].item())
            res = float(cache["residual_norm"].item())
            solver_exit.append(exit_flag)
            residual.append(res)
            unbounded_flag = not torch.isfinite(x_next).all().item()
            unbounded.append(unbounded_flag)
            target_ok = np.isfinite(horizon_targets[0].detach().cpu().numpy()).all()
            keep_flags.append(
                bool(
                    (exit_flag == 0)
                    and (res <= cfg_dyn.residual_threshold)
                    and not unbounded_flag
                    and target_ok
                    and np.isfinite(tip_list[-1]).all()
                )
            )

            if heartbeat_steps and (t % heartbeat_steps == 0):
                print(
                    f"heartbeat mode={mode_name} step={t} elapsed_sec={time.time() - start_time:.2f} "
                    f"solver_exit={exit_flag} residual={res:.3e} unbounded={bool(unbounded_flag)} "
                    f"kept_steps={kept_steps}",
                    flush=True,
                )

            if exit_flag != 0 or res > cfg_dyn.residual_threshold:
                t = n_steps
                break
            if unbounded_flag:
                t = n_steps
                break
            kept_steps += 1
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
        timeout_reason,
        exception_type,
        exception_msg,
        exception_tb,
        np.asarray(keep_flags, dtype=bool),
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


def _first_fail_details(tip, target, solver_exit, residual, unbounded, residual_threshold):
    n = min(tip.shape[0], target.shape[0], solver_exit.shape[0], residual.shape[0], unbounded.shape[0])
    if n == 0:
        return "no_steps", None, None, None, None
    for idx in range(n):
        if unbounded[idx]:
            return "unbounded", idx, int(solver_exit[idx]), float(residual[idx]), True
        if not np.isfinite(tip[idx]).all():
            return "nonfinite_tip", idx, int(solver_exit[idx]), float(residual[idx]), bool(unbounded[idx])
        if not np.isfinite(target[idx]).all():
            return "target_nonfinite", idx, int(solver_exit[idx]), float(residual[idx]), bool(unbounded[idx])
        if solver_exit[idx] != 0:
            return "solver_exit", idx, int(solver_exit[idx]), float(residual[idx]), bool(unbounded[idx])
        if residual[idx] > residual_threshold:
            return "residual", idx, int(solver_exit[idx]), float(residual[idx]), bool(unbounded[idx])
    return "none", None, None, None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=str, default=None)
    parser.add_argument("--dataset", type=str, default="data/dyn_fk_ramp_circle1_hold1.npz")
    parser.add_argument("--start-idx", type=int, default=0)
    parser.add_argument("--end-idx", type=int, default=80)
    parser.add_argument("--window-len", type=int, default=None)
    parser.add_argument("--preset", choices=["acceptance_fast"], default=None)
    parser.add_argument("--init-mode", choices=["dataset_x0", "init_from_ramp", "zero"], default="dataset_x0")
    parser.add_argument("--ramp-len", type=int, default=None)
    parser.add_argument("--horizon", type=int, default=None)
    parser.add_argument("--preroll-steps", type=int, default=None)
    parser.add_argument("--mpc-steps", type=int, default=None)
    parser.add_argument("--max-wall-sec", type=float, default=None)
    parser.add_argument("--heartbeat-steps", type=int, default=None)
    parser.add_argument("--heartbeat-every", type=int, default=None)
    parser.add_argument("--plot", type=int, default=None)
    parser.add_argument("--umax", type=float, default=None)
    parser.add_argument("--d-umax", type=float, default=None)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-prefix", type=str, default=None)
    parser.add_argument("--report-path", type=str, default=None)
    args = parser.parse_args()

    if args.preset == "acceptance_fast":
        if args.window_len is None:
            args.window_len = 80
        if args.preroll_steps is None:
            args.preroll_steps = 10
        if args.horizon is None:
            args.horizon = 5
        if args.mpc_steps is None:
            args.mpc_steps = 5
        if args.max_wall_sec is None:
            args.max_wall_sec = 60.0
        if args.plot is None:
            args.plot = 0
        if args.heartbeat_steps is None:
            args.heartbeat_steps = 5
    if args.horizon is None:
        args.horizon = 15
    if args.preroll_steps is None:
        args.preroll_steps = 0
    if args.mpc_steps is None:
        args.mpc_steps = 1
    if args.max_wall_sec is None:
        args.max_wall_sec = 60.0
    if args.heartbeat_steps is None:
        args.heartbeat_steps = 10
    if args.plot is None:
        args.plot = 0
    if args.heartbeat_every is not None:
        args.heartbeat_steps = int(args.heartbeat_every)

    plot_enabled = bool(args.plot)
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

    timeout_any = False
    timeout_reason = "none"
    timeout_mode = None
    timeout_step = None
    target_seq = np.zeros((0, 3), dtype=np.float64)
    start_idx = args.start_idx
    end_idx = args.end_idx
    requested_start_idx = args.start_idx
    requested_end_idx = args.end_idx
    dt = 0.0
    li_val = 0.0
    umax_safe = 0.0
    d_umax_safe = 0.0
    preroll_x = np.zeros((0, 39), dtype=np.float64)
    preroll_tip = np.zeros((0, 3), dtype=np.float64)
    preroll_solver = np.zeros((0,), dtype=np.int32)
    preroll_residual = np.zeros((0,), dtype=np.float64)
    preroll_unbounded = np.zeros((0,), dtype=bool)
    preroll_timeout = False
    preroll_timeout_step = None
    preroll_u = np.zeros((0, 3), dtype=np.float64)
    prev_u_source = "zero"
    base_start_idx = args.start_idx
    window_len_explicit = args.window_len is not None
    window = 0
    ramp_len = 0
    modes = {}
    metrics = {}
    cfg_max_iter = 1

    dataset_name = os.path.basename(args.dataset)
    n_total = 0
    ext_load_sec = 0.0
    dataset_load_sec = 0.0
    model_load_sec = 0.0
    rollout_elapsed_sec = 0.0

    try:
        if args.model is None:
            model_dir = os.path.join("output_data", "models")
            if not os.path.isdir(model_dir):
                raise SystemExit("--model is required (no models found)")
            candidates = [os.path.join(model_dir, name) for name in os.listdir(model_dir) if name.endswith(".pt")]
            if not candidates:
                raise SystemExit("--model is required (no .pt models found)")
            args.model = max(candidates, key=os.path.getmtime)
            print(f"Using latest model: {args.model}")
    
        t_total_start = time.time()
        print("Loading model...", flush=True)
        t_model_start = time.time()
        model, model_meta = _load_model(args.model)
        model_load_sec = time.time() - t_model_start
        print("Model loaded.", flush=True)
        include_tip = bool(model_meta.get("include_tip", False))
        include_prev_u = bool(model_meta.get("include_prev_u", False))
    
        dataset_path = args.dataset
        if not os.path.exists(dataset_path):
            candidate = os.path.join("data", dataset_path)
            if os.path.exists(candidate):
                dataset_path = candidate
        print("Loading dataset...", flush=True)
        t_data_start = time.time()
        targets, currents, dt, li_val, segments = _load_dataset(dataset_path)
        dataset_load_sec = time.time() - t_data_start
        print("Dataset loaded.", flush=True)
        if currents is None:
            raise SystemExit("dataset currents required")
        dataset_name = os.path.basename(dataset_path)
        n_total = targets.shape[0]
        dataset_payload = np.load(dataset_path)

        ramp_len = _infer_ramp_len(segments, args.ramp_len) if args.init_mode == "init_from_ramp" else 0

        base_start_idx = args.start_idx
        if args.window_len is not None:
            args.end_idx = base_start_idx + int(args.window_len)
            requested_end_idx = args.end_idx
        if not window_len_explicit and args.end_idx == base_start_idx:
            raise ValueError("window_len is 0 without explicit --window-len 0")
        if args.init_mode == "init_from_ramp" and ramp_len > 0:
            base_start_idx = max(base_start_idx, ramp_len)

        preroll_steps = max(int(args.preroll_steps), 0)
        window = max(args.end_idx - base_start_idx, 0)
        if window == 0 and not window_len_explicit:
            raise ValueError("window_len is 0 without explicit --window-len 0")
        start_idx = base_start_idx + preroll_steps
        end_idx = start_idx + window
        target_seq = _prepare_targets(targets, start_idx, window)

        print("Loading dynamics extension...", flush=True)
        t_ext_start = time.time()
        ext = load_dyn_ext()
        ext_load_sec = time.time() - t_ext_start
        print("Dynamics extension loaded.", flush=True)
        cfg_dyn = build_cfg(ext)
        x0_default, _ = build_state()

        safe_umax, safe_d_umax = _load_safe_bounds_optional()
        dataset_umax = dataset_payload.get("umax")
        dataset_d_umax = dataset_payload.get("d_umax")
        if args.umax is not None:
            umax_safe = float(args.umax)
        elif safe_umax is not None:
            umax_safe = float(safe_umax)
        elif dataset_umax is not None:
            umax_safe = float(dataset_umax)
        else:
            umax_safe = 0.1
        if args.d_umax is not None:
            d_umax_safe = float(args.d_umax)
        elif safe_d_umax is not None:
            d_umax_safe = float(safe_d_umax)
        elif dataset_d_umax is not None:
            d_umax_safe = float(dataset_d_umax)
        else:
            d_umax_safe = 0.01
        if args.preset and args.umax is None and umax_safe <= 0.0:
            raise ValueError("umax must be >0 for presets unless explicitly set to 0")
        if args.preset and args.d_umax is None and d_umax_safe <= 0.0:
            raise ValueError("d_umax must be >0 for presets unless explicitly set to 0")

        li = torch.tensor([li_val], dtype=torch.float64)
    
        x0 = x0_default.clone()
        prefix_timeout = False
        prefix_timeout_step = None
        if args.init_mode == "init_from_ramp" and ramp_len > 0:
            ramp_u = _clamp_currents(currents[:ramp_len], umax_safe, d_umax_safe)
            x0, prefix_timeout, prefix_timeout_step = _rollout_prefix(
                x0, ramp_u, li, cfg_dyn, mode_name="prefix_ramp"
            )
            if base_start_idx > ramp_len and not prefix_timeout:
                prefix_u = _clamp_currents(currents[ramp_len:base_start_idx], umax_safe, d_umax_safe)
                x0, prefix_timeout, prefix_timeout_step = _rollout_prefix(
                    x0, prefix_u, li, cfg_dyn, mode_name="prefix_base"
                )
        elif args.init_mode == "dataset_x0" and base_start_idx > 0:
            prefix_u = _clamp_currents(currents[:base_start_idx], umax_safe, d_umax_safe)
            x0, prefix_timeout, prefix_timeout_step = _rollout_prefix(
                x0, prefix_u, li, cfg_dyn, mode_name="prefix_dataset"
            )
    
        preroll_start_idx = base_start_idx
        preroll_end_idx = base_start_idx + preroll_steps
        preroll_u = _clamp_currents(currents[preroll_start_idx:preroll_end_idx], umax_safe, d_umax_safe)
        preroll_timeout = False
        preroll_timeout_step = None
        if preroll_steps > 0 and preroll_u.shape[0] > 0:
            (
                preroll_x,
                preroll_tip,
                preroll_solver,
                preroll_residual,
                preroll_unbounded,
                x0,
                preroll_timeout,
                preroll_timeout_step,
            ) = _rollout_open_loop(
                x0,
                preroll_u,
                li,
                cfg_dyn,
                mode_name="preroll",
            )
        else:
            preroll_x = np.zeros((0, 39), dtype=np.float64)
            preroll_tip = np.zeros((0, 3), dtype=np.float64)
            preroll_solver = np.zeros((0,), dtype=np.int32)
            preroll_residual = np.zeros((0,), dtype=np.float64)
            preroll_unbounded = np.zeros((0,), dtype=bool)
        if preroll_timeout:
            prefix_timeout = True
            prefix_timeout_step = preroll_timeout_step
    
        window = max(args.end_idx - base_start_idx, 0)
        start_idx = base_start_idx + preroll_steps
        end_idx = start_idx + window
        target_seq = _prepare_targets(targets, start_idx, window)
        n_total = targets.shape[0]
        expert_u = _clamp_currents(currents[start_idx : start_idx + window], umax_safe, d_umax_safe)
        if preroll_u.shape[0]:
            prev_u_init = preroll_u[-1].copy()
            prev_u_source = "preroll_u_last"
        else:
            prev_u_init = expert_u[0] if expert_u.shape[0] else np.zeros((3,), dtype=np.float64)
            prev_u_source = "dataset_u0" if expert_u.shape[0] else "zero"
        print(f"prev_u_source={prev_u_source}", flush=True)
        
        def _empty_mode(reason="max_wall_sec"):
            return (
                np.zeros((0, 39), dtype=np.float64),
                np.zeros((0, 3), dtype=np.float64),
                np.zeros((0, 3), dtype=np.float64),
                np.zeros((0, 3), dtype=np.float64),
                np.zeros((0,), dtype=np.int32),
                np.zeros((0,), dtype=np.float64),
                np.zeros((0,), dtype=bool),
                True,
                0,
                reason,
            )
    
        results = {}
        if window < 2:
            results["mpc_only"] = _empty_mode("short_window")
        if args.mpc_steps is not None and int(args.mpc_steps) < 1:
            results["mpc_only"] = _empty_mode("invalid_mpc_steps")
        if prefix_timeout:
            results["policy_only"] = _empty_mode("step_timeout")
            results["mpc_only"] = _empty_mode("step_timeout")
            results["policy_plus_mpc"] = _empty_mode("step_timeout")
            policy_runtime = 0.0
            mpc_runtime = 0.0
            pol_mpc_runtime = 0.0
        
        cfg_max_iter = 1
        rollout_elapsed_holder = [0.0]
        t_rollout_start = time.time()
        eval_start = t_rollout_start
        if not prefix_timeout:
            t0 = eval_start
            pol = _rollout_policy_only_with_walltime(
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
            t_rollout_start,
            args.heartbeat_steps,
            "policy_only",
            rollout_elapsed_holder,
        )
            results["policy_only"] = pol
            policy_runtime = time.time() - t0
    
            if pol[7] or (time.time() - eval_start >= args.max_wall_sec):
                results["mpc_only"] = _empty_mode()
                results["policy_plus_mpc"] = _empty_mode()
                mpc_runtime = 0.0
                pol_mpc_runtime = 0.0
            else:
                cfg = ILQRConfig(
            horizon=args.horizon,
            max_iter=cfg_max_iter,
                    w_tip=1.0,
                    w_tip_terminal=10.0,
                    w_u=1e-2,
                    w_du=1e-1,
                    max_du=d_umax_safe,
                    max_u=umax_safe,
                    linearization_backend="v1_3",
                    linearization_dense=True,
                )
    
                def _zeros_warm(_x, _targets):
                    return torch.zeros((cfg.horizon, 1, 3), dtype=torch.float64)
    
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
    
                t0 = time.time()
                mpc_only = _mpc_loop(
                    x0,
                    torch.tensor(target_seq, dtype=torch.float64),
                    li,
                    cfg_dyn,
                    cfg,
                    _zeros_warm,
                    args.mpc_steps,
            args.max_wall_sec,
            t_rollout_start,
            args.heartbeat_steps,
            "mpc_only",
            rollout_elapsed_holder,
        )
                mpc_runtime = time.time() - t0
    
                t0 = time.time()
                pol_mpc = _mpc_loop(
                    x0,
                    torch.tensor(target_seq, dtype=torch.float64),
                    li,
                    cfg_dyn,
                    cfg,
                    _policy_warm,
                    args.mpc_steps,
            args.max_wall_sec,
            t_rollout_start,
            args.heartbeat_steps,
            "policy_plus_mpc",
            rollout_elapsed_holder,
        )
                pol_mpc_runtime = time.time() - t0
                results["mpc_only"] = mpc_only
                results["policy_plus_mpc"] = pol_mpc
    
        pol = results["policy_only"]
        mpc_only = results["mpc_only"]
        pol_mpc = results["policy_plus_mpc"]
    
        modes = {
            "policy_only": {
                "x": pol[0],
                "tip": pol[1],
                "u": pol[2],
                "u_raw": pol[3],
                "solver_exit": pol[4],
                "residual": pol[5],
                "unbounded": pol[6],
                "timeout": pol[7],
                "timeout_step": pol[8],
                "timeout_reason": pol[9],
                "runtime": policy_runtime,
            },
            "mpc_only": {
                "x": mpc_only[0],
                "tip": mpc_only[1],
                "u": mpc_only[2],
                "solver_exit": mpc_only[3],
                "residual": mpc_only[4],
                "unbounded": mpc_only[5],
                "timeout": mpc_only[6],
                "timeout_step": mpc_only[7],
                "timeout_reason": mpc_only[8],
                "exception_type": mpc_only[9] if len(mpc_only) > 9 else "",
                "exception_msg": mpc_only[10] if len(mpc_only) > 10 else "",
                "exception_tb": mpc_only[11] if len(mpc_only) > 11 else "",
                "keep_override": mpc_only[12] if len(mpc_only) > 12 else None,
                "runtime": mpc_runtime,
            },
            "policy_plus_mpc": {
                "x": pol_mpc[0],
                "tip": pol_mpc[1],
                "u": pol_mpc[2],
                "solver_exit": pol_mpc[3],
                "residual": pol_mpc[4],
                "unbounded": pol_mpc[5],
                "timeout": pol_mpc[6],
                "timeout_step": pol_mpc[7],
                "timeout_reason": pol_mpc[8],
                "exception_type": pol_mpc[9] if len(pol_mpc) > 9 else "",
                "exception_msg": pol_mpc[10] if len(pol_mpc) > 10 else "",
                "exception_tb": pol_mpc[11] if len(pol_mpc) > 11 else "",
                "runtime": pol_mpc_runtime,
            },
        }
    
        metrics = {}
        for name, payload in modes.items():
            if payload.get("keep_override") is not None:
                keep = payload["keep_override"]
            else:
                keep = _keep_mask(
                    payload["tip"],
                    target_seq,
                    payload["solver_exit"],
                    payload["residual"],
                    payload["unbounded"],
                    cfg_dyn.residual_threshold,
                )
            n_err = min(payload["tip"].shape[0], target_seq.shape[0])
            err = (
                np.linalg.norm(payload["tip"][:n_err] - target_seq[:n_err], axis=1)
                if n_err
                else np.zeros((0,))
            )
            err_finite = np.isfinite(err)
            if err.shape[0] != keep.shape[0]:
                keep = keep[: err.shape[0]]
            keep = keep & err_finite
            mean_err, max_err, _ = _error_metrics(payload["tip"], target_seq, keep)
            metrics[name] = {
                "keep": keep,
                "mean_err": mean_err,
                "max_err": max_err,
                "err": err[keep] if keep.any() else np.zeros((0,)),
                "attempted": int(target_seq.shape[0]),
                "kept": int(keep.sum()),
                "dropped": int(keep.shape[0] - keep.sum()),
                "runtime": payload["runtime"],
                "timeout": bool(payload["timeout"]),
                "timeout_step": payload["timeout_step"],
                "timeout_reason": payload["timeout_reason"],
            }
    
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_prefix = args.output_prefix or f"policy_warmstart_{stamp}"
        out_path = os.path.join("output_data", f"{out_prefix}.npz")
    
        timeout_any = any(metrics[name]["timeout"] for name in metrics)
        if timeout_any:
            plot_enabled = False

        _atomic_save_npz(
            out_path,
            target_xyz=target_seq,
            start_idx=start_idx,
            end_idx=end_idx,
            dt=np.asarray(dt, dtype=np.float64),
            Li=np.asarray(li_val, dtype=np.float64),
            umax=np.asarray(umax_safe, dtype=np.float64),
            d_umax=np.asarray(d_umax_safe, dtype=np.float64),
            preroll_steps=preroll_steps,
            preroll_start_idx=preroll_start_idx,
            preroll_end_idx=preroll_end_idx,
            x_preroll=preroll_x,
            tip_preroll=preroll_tip,
            u_preroll=preroll_u,
            solver_exit_preroll=preroll_solver,
            residual_preroll=preroll_residual,
            unbounded_preroll=preroll_unbounded,
            timeout_preroll=preroll_timeout,
            timeout_step_preroll=preroll_timeout_step,
            max_wall_sec=args.max_wall_sec,
            mpc_steps=args.mpc_steps,
            **{f"tip_{name}": modes[name]["tip"] for name in modes},
            **{f"u_{name}": modes[name]["u"] for name in modes},
            **{f"solver_exit_{name}": modes[name]["solver_exit"] for name in modes},
            **{f"residual_{name}": modes[name]["residual"] for name in modes},
            **{f"unbounded_{name}": modes[name]["unbounded"] for name in modes},
            **{f"keep_{name}": metrics[name]["keep"] for name in metrics},
            **{f"error_{name}": metrics[name]["err"] for name in metrics},
            **{f"timeout_{name}": metrics[name]["timeout"] for name in metrics},
            **{f"timeout_step_{name}": metrics[name]["timeout_step"] for name in metrics},
        )
    
        if plot_enabled:
            try:
                import matplotlib.pyplot as plt
            except ImportError:
                plt = None
            if plt is not None:
                fig, ax = plt.subplots(figsize=(8, 6))
                for name in ["policy_only", "mpc_only", "policy_plus_mpc"]:
                    err = metrics[name]["err"]
                    ax.plot(err, label=name)
                ax.set_title("Kept error per mode")
                ax.set_xlabel("time")
                ax.set_ylabel("error")
                ax.legend()
                fig.tight_layout()
                out_plot = os.path.join("output_data", f"{out_prefix}_error")
                fig.savefig(out_plot + ".png", dpi=300)
                fig.savefig(out_plot + ".pdf")

            plot_utils.save_tip_trajectory_plot(
                modes["policy_only"]["tip"],
                target_seq[: modes["policy_only"]["tip"].shape[0]],
                title="Policy Only",
                out_path=os.path.join("output_data", f"{out_prefix}_policy_only"),
            )
            plot_utils.save_tip_trajectory_plot(
                modes["mpc_only"]["tip"],
                target_seq[: modes["mpc_only"]["tip"].shape[0]],
                title="MPC Only",
                out_path=os.path.join("output_data", f"{out_prefix}_mpc_only"),
            )
            plot_utils.save_tip_trajectory_plot(
                modes["policy_plus_mpc"]["tip"],
                target_seq[: modes["policy_plus_mpc"]["tip"].shape[0]],
                title="Policy + MPC",
                out_path=os.path.join("output_data", f"{out_prefix}_policy_plus_mpc"),
            )
            plot_prefix = os.path.join("output_data", out_prefix)
        else:
            plot_prefix = None
    
        report_path = args.report_path or os.path.join(
            "docs",
            "control",
            "run_reports",
            f"policy_warmstart_mpc_{stamp}.md",
        )
    
        rollout_elapsed_sec = rollout_elapsed_holder[0]
        total_elapsed_sec = time.time() - t_total_start
        effective_cfg = {
            "preset": args.preset or "none",
            "start_idx": start_idx,
            "end_idx": end_idx,
            "window_len": window,
            "horizon": args.horizon,
            "mpc_steps": args.mpc_steps,
            "max_iter": cfg_max_iter,
            "max_wall_sec": args.max_wall_sec,
            "preroll_steps": preroll_steps,
            "init_mode": args.init_mode,
            "plot": bool(plot_enabled),
            "seed": args.seed,
            "heartbeat_steps": args.heartbeat_steps,
        }
        print("Effective config:")
        print(json.dumps(effective_cfg, indent=2))
    
        timeout_any = any(metrics[name]["timeout"] for name in metrics)
        eval_go = bool(metrics["policy_plus_mpc"]["kept"] > 0) and not timeout_any
        mpc_fail_reason = ""
        mpc_fail_step = None
        mpc_fail_solver = None
        mpc_fail_residual = None
        mpc_fail_unbounded = None
        if metrics["mpc_only"]["kept"] == 0:
            (
                mpc_fail_reason,
                mpc_fail_step,
                mpc_fail_solver,
                mpc_fail_residual,
                mpc_fail_unbounded,
            ) = _first_fail_details(
                modes["mpc_only"]["tip"],
                target_seq,
                modes["mpc_only"]["solver_exit"],
                modes["mpc_only"]["residual"],
                modes["mpc_only"]["unbounded"],
                cfg_dyn.residual_threshold,
            )

        _atomic_write_report(
            report_path=report_path,
            command=" ".join(sys.argv),
            dataset_name=dataset_name,
            start_idx=start_idx,
            end_idx=end_idx,
            n_total=n_total,
            li_mm=li_val,
            dt=dt,
            umax=umax_safe,
            d_umax=d_umax_safe,
            mean_error=metrics["policy_plus_mpc"]["mean_err"],
            max_error=metrics["policy_plus_mpc"]["max_err"],
            artifact_npz=out_path,
            plot_prefix=plot_prefix,
            init_from_ramp=args.init_mode == "init_from_ramp",
            ramp_len=ramp_len,
            max_wall_hit=bool(timeout_any),
            extra_entries={
                "requested_start_idx": requested_start_idx,
                "requested_end_idx": requested_end_idx,
                "init_mode": args.init_mode,
                "prev_u_source": prev_u_source,
                "base_start_idx": base_start_idx,
                "preroll_steps": preroll_steps,
                "preroll_start_idx": preroll_start_idx,
                "preroll_end_idx": preroll_end_idx,
                "max_wall_sec": args.max_wall_sec,
                "mpc_steps": args.mpc_steps,
                "ext_load_sec": ext_load_sec,
                "dataset_load_sec": dataset_load_sec,
                "model_load_sec": model_load_sec,
                "rollout_elapsed_sec": rollout_elapsed_sec,
                "total_elapsed_sec": total_elapsed_sec,
                "effective_config": effective_cfg,
                "policy_only_mean": metrics["policy_only"]["mean_err"],
                "policy_only_max": metrics["policy_only"]["max_err"],
                "mpc_only_mean": metrics["mpc_only"]["mean_err"],
                "mpc_only_max": metrics["mpc_only"]["max_err"],
                "policy_plus_mpc_mean": metrics["policy_plus_mpc"]["mean_err"],
                "policy_plus_mpc_max": metrics["policy_plus_mpc"]["max_err"],
                "policy_only_kept": metrics["policy_only"]["kept"],
                "mpc_only_kept": metrics["mpc_only"]["kept"],
                "policy_plus_mpc_kept": metrics["policy_plus_mpc"]["kept"],
                "policy_only_runtime_sec": metrics["policy_only"]["runtime"],
                "mpc_only_runtime_sec": metrics["mpc_only"]["runtime"],
                "policy_plus_mpc_runtime_sec": metrics["policy_plus_mpc"]["runtime"],
                "policy_only_timeout": metrics["policy_only"]["timeout"],
                "mpc_only_timeout": metrics["mpc_only"]["timeout"],
                "policy_plus_mpc_timeout": metrics["policy_plus_mpc"]["timeout"],
                "policy_only_timeout_step": metrics["policy_only"]["timeout_step"],
                "mpc_only_timeout_step": metrics["mpc_only"]["timeout_step"],
                "policy_plus_mpc_timeout_step": metrics["policy_plus_mpc"]["timeout_step"],
                "policy_only_timeout_reason": metrics["policy_only"]["timeout_reason"],
                "mpc_only_timeout_reason": metrics["mpc_only"]["timeout_reason"],
                "policy_plus_mpc_timeout_reason": metrics["policy_plus_mpc"]["timeout_reason"],
                "mpc_only_exception_type": modes["mpc_only"].get("exception_type", ""),
                "mpc_only_exception_msg": modes["mpc_only"].get("exception_msg", ""),
                "mpc_only_exception_tb": modes["mpc_only"].get("exception_tb", ""),
                "failure_reason_mpc_only": mpc_fail_reason,
                "first_fail_step_mpc_only": mpc_fail_step if mpc_fail_step is not None else "n/a",
                "first_fail_solver_exit_mpc_only": mpc_fail_solver if mpc_fail_solver is not None else "n/a",
                "first_fail_residual_mpc_only": mpc_fail_residual if mpc_fail_residual is not None else "n/a",
                "first_fail_unbounded_mpc_only": mpc_fail_unbounded if mpc_fail_unbounded is not None else "n/a",
                "eval_go": bool(eval_go),
            },
        )
    
        print(f"Saved {out_path}")
        print(f"Report {report_path}")
        try:
            np.load(out_path, allow_pickle=False)
            print("NPZ load check: allow_pickle=False OK", flush=True)
        except Exception as exc:
            print(f"NPZ load check failed: {exc}", flush=True)
        if timeout_any:
            print("Timed out (max_wall_sec reached); exiting with code 2", flush=True)
            raise SystemExit(2)
    except Exception as exc:
        print(f"Exception in eval: {exc}", flush=True)
        print(traceback.format_exc(), flush=True)
        timeout_any = False
        timeout_reason = "exception"
        timeout_mode = "global"
        timeout_step = None
        modes = {
            "policy_only": {
                "tip": np.zeros((0, 3), dtype=np.float64),
                "u": np.zeros((0, 3), dtype=np.float64),
                "u_raw": np.zeros((0, 3), dtype=np.float64),
                "solver_exit": np.zeros((0,), dtype=np.int32),
                "residual": np.zeros((0,), dtype=np.float64),
                "unbounded": np.zeros((0,), dtype=bool),
                "timeout": True,
                "timeout_step": 0,
                "timeout_reason": timeout_reason,
                "runtime": 0.0,
            },
            "mpc_only": {
                "tip": np.zeros((0, 3), dtype=np.float64),
                "u": np.zeros((0, 3), dtype=np.float64),
                "solver_exit": np.zeros((0,), dtype=np.int32),
                "residual": np.zeros((0,), dtype=np.float64),
                "unbounded": np.zeros((0,), dtype=bool),
                "timeout": True,
                "timeout_step": 0,
                "timeout_reason": timeout_reason,
                "runtime": 0.0,
            },
            "policy_plus_mpc": {
                "tip": np.zeros((0, 3), dtype=np.float64),
                "u": np.zeros((0, 3), dtype=np.float64),
                "solver_exit": np.zeros((0,), dtype=np.int32),
                "residual": np.zeros((0,), dtype=np.float64),
                "unbounded": np.zeros((0,), dtype=bool),
                "timeout": True,
                "timeout_step": 0,
                "timeout_reason": timeout_reason,
                "runtime": 0.0,
            },
        }
        metrics = {
            name: {
                "keep": np.zeros((0,), dtype=bool),
                "mean_err": float("nan"),
                "max_err": float("nan"),
                "err": np.zeros((0,)),
                "attempted": 0,
                "kept": 0,
                "dropped": 0,
                "runtime": 0.0,
                "timeout": True,
                "timeout_step": 0,
                "timeout_reason": timeout_reason,
            }
            for name in modes
        }
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_prefix = args.output_prefix or f"policy_warmstart_{stamp}"
        out_path = os.path.join("output_data", f"{out_prefix}.npz")
        report_path = args.report_path or os.path.join(
            "docs",
            "control",
            "run_reports",
            f"policy_warmstart_mpc_{stamp}.md",
        )
        _atomic_save_npz(
            out_path,
            target_xyz=target_seq,
            start_idx=start_idx,
            end_idx=end_idx,
            dt=np.asarray(dt, dtype=np.float64),
            Li=np.asarray(li_val, dtype=np.float64),
            umax=np.asarray(umax_safe, dtype=np.float64),
            d_umax=np.asarray(d_umax_safe, dtype=np.float64),
            preroll_steps=int(args.preroll_steps),
            preroll_start_idx=base_start_idx,
            preroll_end_idx=base_start_idx + int(args.preroll_steps),
            x_preroll=preroll_x,
            tip_preroll=preroll_tip,
            u_preroll=preroll_u,
            solver_exit_preroll=preroll_solver,
            residual_preroll=preroll_residual,
            unbounded_preroll=preroll_unbounded,
            timeout_preroll=preroll_timeout,
            timeout_step_preroll=preroll_timeout_step,
            max_wall_sec=args.max_wall_sec,
            mpc_steps=args.mpc_steps,
            **{f"tip_{name}": modes[name]["tip"] for name in modes},
            **{f"u_{name}": modes[name]["u"] for name in modes},
            **{f"solver_exit_{name}": modes[name]["solver_exit"] for name in modes},
            **{f"residual_{name}": modes[name]["residual"] for name in modes},
            **{f"unbounded_{name}": modes[name]["unbounded"] for name in modes},
            **{f"keep_{name}": metrics[name]["keep"] for name in metrics},
            **{f"error_{name}": metrics[name]["err"] for name in metrics},
            **{f"timeout_{name}": metrics[name]["timeout"] for name in metrics},
            **{f"timeout_step_{name}": metrics[name]["timeout_step"] for name in metrics},
        )
        effective_cfg = {
            "preset": args.preset or "none",
            "start_idx": start_idx,
            "end_idx": end_idx,
            "window_len": window,
            "horizon": args.horizon,
            "mpc_steps": args.mpc_steps,
            "max_iter": cfg_max_iter,
            "max_wall_sec": args.max_wall_sec,
            "preroll_steps": int(args.preroll_steps),
            "init_mode": args.init_mode,
            "plot": bool(plot_enabled),
            "seed": args.seed,
            "heartbeat_steps": args.heartbeat_steps,
        }
        total_elapsed_sec = time.time() - t_total_start if "t_total_start" in locals() else 0.0
        _atomic_write_report(
            report_path=report_path,
            command=" ".join(sys.argv),
            dataset_name=os.path.basename(dataset_path) if "dataset_path" in locals() else os.path.basename(args.dataset),
            start_idx=start_idx,
            end_idx=end_idx,
            n_total=n_total if "n_total" in locals() else 0,
            li_mm=li_val,
            dt=dt,
            umax=umax_safe,
            d_umax=d_umax_safe,
            mean_error=float("nan"),
            max_error=float("nan"),
            artifact_npz=out_path,
            plot_prefix=None,
            init_from_ramp=args.init_mode == "init_from_ramp",
            ramp_len=ramp_len,
            max_wall_hit=False,
            extra_entries={
                "requested_start_idx": requested_start_idx,
                "requested_end_idx": requested_end_idx,
                "timeout_reason": timeout_reason,
                "timeout_mode": timeout_mode,
                "timeout_step": timeout_step,
                "ext_load_sec": ext_load_sec,
                "dataset_load_sec": dataset_load_sec,
                "model_load_sec": model_load_sec,
                "rollout_elapsed_sec": rollout_elapsed_sec,
                "total_elapsed_sec": total_elapsed_sec,
                "effective_config": effective_cfg,
                "eval_go": False,
            },
        )
        print(f"Saved {out_path}")
        print(f"Report {report_path}")
        try:
            np.load(out_path, allow_pickle=False)
            print("NPZ load check: allow_pickle=False OK", flush=True)
        except Exception as exc:
            print(f"NPZ load check failed: {exc}", flush=True)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
