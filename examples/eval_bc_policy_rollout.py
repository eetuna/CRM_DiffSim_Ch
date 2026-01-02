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


def _prepare_targets(targets, start_idx, n_steps):
    end_idx = start_idx + n_steps
    clipped = targets[start_idx:end_idx]
    if clipped.shape[0] < n_steps:
        pad = np.repeat(clipped[-1:], n_steps - clipped.shape[0], axis=0)
        clipped = np.concatenate([clipped, pad], axis=0)
    return clipped


def _rollout_expert(x0, u_seq, targets, li, cfg_dyn):
    x_list = []
    tip_list = []
    solver_exit = []
    residual_norm = []
    unbounded = []

    x_t = x0
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
        for t in range(u_seq.shape[0]):
            x_list.append(x_t.detach().cpu().numpy()[0])
            tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())
            u_t = torch.tensor(u_seq[t : t + 1], dtype=torch.float64).view(1, 1, 3)
            x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
            exit_flag = int(cache["solver_exit"].item())
            residual = float(cache["residual_norm"].item())
            solver_exit.append(exit_flag)
            residual_norm.append(residual)
            message = stdout_buf.getvalue() + stderr_buf.getvalue()
            unbounded_flag = "Unbounded" in message
            if not torch.isfinite(x_tp1).all().item():
                unbounded_flag = True
            unbounded.append(unbounded_flag)

            if exit_flag != 0 or residual > cfg_dyn.residual_threshold:
                break
            if unbounded_flag or not torch.isfinite(x_tp1).all().item():
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
        errors,
    )


def _rollout_policy(
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
    debug_first_step,
):
    x_list = []
    tip_list = []
    u_list = []
    u_raw_list = []
    solver_exit = []
    residual_norm = []
    unbounded = []

    x_t = x0
    u_prev = prev_u_init.copy()
    debug_payload = None
    stdout_buf = io.StringIO()
    stderr_buf = io.StringIO()
    with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
        for t in range(targets.shape[0]):
            x_list.append(x_t.detach().cpu().numpy()[0])
            tip_list.append(x_t[0, TIP_IDX].detach().cpu().numpy())

            if include_tip:
                state = x_t[0, TIP_IDX].detach().cpu().float()
            else:
                state = x_t[0].detach().cpu().float()
            target = torch.tensor(targets[t], dtype=torch.float32)
            inputs = [state, target]
            if include_prev_u:
                inputs.append(torch.tensor(u_prev, dtype=torch.float32))
            x_in = torch.cat(inputs, dim=0).unsqueeze(0)

            with torch.no_grad():
                u_pred = model(x_in)[0].cpu().numpy()
            u_raw_list.append(u_pred.copy())
            u_clamped = np.clip(u_pred, -umax, umax)
            delta = u_clamped - u_prev
            delta = np.clip(delta, -d_umax, d_umax)
            u_clamped = u_prev + delta
            u_prev = u_clamped.copy()
            u_list.append(u_clamped)

            u_t = torch.tensor(u_clamped, dtype=torch.float64).view(1, 1, 3)
            x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
            exit_flag = int(cache["solver_exit"].item())
            residual = float(cache["residual_norm"].item())
            solver_exit.append(exit_flag)
            residual_norm.append(residual)

            message = stdout_buf.getvalue() + stderr_buf.getvalue()
            unbounded_flag = "Unbounded" in message
            if not torch.isfinite(x_tp1).all().item():
                unbounded_flag = True
            unbounded.append(unbounded_flag)

            if debug_first_step and t == 0:
                debug_payload = {
                    "x0_finite": bool(torch.isfinite(x_t).all().item()),
                    "x0_tip": x_t[0, TIP_IDX].detach().cpu().numpy(),
                    "target0": targets[0],
                    "u_pred_raw": u_pred,
                    "u_pred_clamped": u_clamped,
                    "u_norm": float(np.linalg.norm(u_pred)),
                    "du_norm": float(np.linalg.norm(delta)),
                    "solver_exit": exit_flag,
                    "residual_norm": residual,
                    "unbounded": bool(unbounded_flag),
                }

            if exit_flag != 0 or residual > cfg_dyn.residual_threshold:
                break
            if unbounded_flag or not torch.isfinite(x_tp1).all().item():
                break
            x_t = x_tp1

    tip_np = np.asarray(tip_list, dtype=np.float64)
    target_np = targets[: tip_np.shape[0]]
    errors = np.linalg.norm(tip_np - target_np, axis=1) if tip_np.size else np.zeros((0,))

    return (
        np.asarray(x_list, dtype=np.float64),
        tip_np,
        np.asarray(u_list, dtype=np.float64),
        np.asarray(u_raw_list, dtype=np.float64),
        np.asarray(solver_exit, dtype=np.int32),
        np.asarray(residual_norm, dtype=np.float64),
        np.asarray(unbounded, dtype=bool),
        errors,
        debug_payload,
    )


def _plot(out_dir, prefix, target_xyz, tip_policy, tip_expert, err_policy, err_expert):
    if tip_policy.size:
        plot_utils.save_tip_trajectory_plot(
            tip_policy,
            target_xyz,
            title="Policy rollout",
            out_path=os.path.join(out_dir, f"{prefix}_tip"),
            planned_xyz=tip_expert,
        )
    if err_policy.size:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            return
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(err_policy, label="policy")
        ax.plot(err_expert, label="expert")
        ax.set_title("Tip error")
        ax.set_xlabel("time")
        ax.set_ylabel("error")
        ax.legend()
        fig.tight_layout()
        out_path = os.path.join(out_dir, f"{prefix}_error")
        fig.savefig(out_path + ".png", dpi=300)
        fig.savefig(out_path + ".pdf")


def _load_model(path):
    payload = torch.load(path, map_location="cpu")
    hidden_sizes = payload["hidden_sizes"]
    input_dim = payload["input_dim"]
    model = MLPPolicy(input_dim, hidden_sizes)
    model.load_state_dict(payload["state_dict"])
    model.eval()
    return model, payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", type=str, default=None)
    parser.add_argument("--circle-dataset", type=str, default="data/dyn_fk_ramp_circle1_hold1.npz")
    parser.add_argument("--lemniscate-dataset", type=str, default="data/dyn_fk_lem1_y40_a10_L94_hold1.npz")
    parser.add_argument("--bc-dataset", type=str, default=None)
    parser.add_argument("--start-idx", type=int, default=0)
    parser.add_argument("--subset-len", type=int, default=60)
    parser.add_argument("--init-mode", choices=["dataset_x0", "init_from_ramp", "zero"], default="dataset_x0")
    parser.add_argument("--ramp-len", type=int, default=None)
    parser.add_argument("--preroll-steps", type=int, default=0)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--output-prefix", type=str, default=None)
    parser.add_argument("--report-prefix", type=str, default=None)
    parser.add_argument("--preset", choices=["circle_easy"], default=None)
    parser.add_argument("--debug-first-step", action="store_true")
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

    umax_safe, d_umax_safe = _load_safe_bounds()

    if args.preset == "circle_easy":
        args.circle_dataset = "data/dyn_fk_ramp_circle1_hold1.npz"
        args.start_idx = 0
        args.subset_len = 200

    bc_counts = None
    if args.bc_dataset:
        data = np.load(args.bc_dataset, allow_pickle=True)
        drop_reasons = data.get("drop_reasons_counts", "{}")
        if isinstance(drop_reasons, np.ndarray):
            drop_reasons = drop_reasons.item()
        bc_counts = {
            "n_steps_total": int(data.get("n_steps_total", 0)),
            "n_steps_kept": int(data.get("n_steps_kept", 0)),
            "n_steps_dropped": int(data.get("n_steps_dropped", 0)),
            "drop_reasons_counts": str(drop_reasons),
        }
        print(
            "dataset_counts: "
            f"total={bc_counts['n_steps_total']} kept={bc_counts['n_steps_kept']} "
            f"dropped={bc_counts['n_steps_dropped']} reasons={bc_counts['drop_reasons_counts']}"
        )

    datasets = [
        ("circle", args.circle_dataset),
        ("lemniscate", args.lemniscate_dataset),
    ]

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for name, path in datasets:
        targets, currents, dt, li_val, segments = _load_dataset(path)
        if currents is None:
            raise SystemExit(f"Dataset {path} missing currents for expert comparison")

        li = torch.tensor([li_val], dtype=torch.float64)
        ramp_len = _infer_ramp_len(segments, args.ramp_len) if args.init_mode == "init_from_ramp" else 0
        base_start_idx = args.start_idx
        if args.init_mode == "init_from_ramp" and ramp_len > 0:
            base_start_idx = max(base_start_idx, ramp_len)

        x0 = x0_default.clone()
        prefix_failed = None
        try:
            if args.init_mode == "init_from_ramp" and ramp_len > 0:
                ramp_u = _clamp_currents(currents[:ramp_len], umax_safe, d_umax_safe)
                x0 = _rollout_prefix(x0, ramp_u, li, cfg_dyn)
                if base_start_idx > ramp_len:
                    prefix_u = _clamp_currents(currents[ramp_len:base_start_idx], umax_safe, d_umax_safe)
                    x0 = _rollout_prefix(x0, prefix_u, li, cfg_dyn)
            elif args.init_mode == "dataset_x0" and base_start_idx > 0:
                prefix_u = _clamp_currents(currents[:base_start_idx], umax_safe, d_umax_safe)
                x0 = _rollout_prefix(x0, prefix_u, li, cfg_dyn)
        except RuntimeError as exc:
            prefix_failed = str(exc)

        preroll_steps = max(int(args.preroll_steps), 0)
        preroll_start_idx = base_start_idx
        preroll_end_idx = base_start_idx + preroll_steps
        preroll_u = _clamp_currents(currents[preroll_start_idx:preroll_end_idx], umax_safe, d_umax_safe)
        if preroll_steps > 0 and preroll_u.shape[0] > 0 and prefix_failed is None:
            try:
                x0 = _rollout_prefix(x0, preroll_u, li, cfg_dyn)
            except RuntimeError as exc:
                prefix_failed = str(exc)

        start_idx = base_start_idx + preroll_steps
        target_seq = _prepare_targets(targets, start_idx, args.subset_len)
        expert_u = _clamp_currents(currents[start_idx : start_idx + args.subset_len], umax_safe, d_umax_safe)
        if preroll_u.shape[0]:
            prev_u_init = preroll_u[-1].copy()
            prev_u_source = "preroll_u_last"
        else:
            prev_u_init = expert_u[0] if expert_u.shape[0] else np.zeros((3,), dtype=np.float64)
            prev_u_source = "dataset_u0" if expert_u.shape[0] else "zero"

        if prefix_failed:
            x_pol = np.zeros((0, 39), dtype=np.float64)
            tip_pol = np.zeros((0, 3), dtype=np.float64)
            u_pol = np.zeros((0, 3), dtype=np.float64)
            u_raw_pol = np.zeros((0, 3), dtype=np.float64)
            solver_pol = np.zeros((0,), dtype=np.int32)
            res_pol = np.zeros((0,), dtype=np.float64)
            unb_pol = np.zeros((0,), dtype=bool)
            tip_exp = np.zeros((0, 3), dtype=np.float64)
            solver_exp = np.zeros((0,), dtype=np.int32)
            res_exp = np.zeros((0,), dtype=np.float64)
            unb_exp = np.zeros((0,), dtype=bool)
            debug_payload = None
        else:
            x_pol, tip_pol, u_pol, u_raw_pol, solver_pol, res_pol, unb_pol, _, debug_payload = _rollout_policy(
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
                args.debug_first_step,
            )
            x_exp, tip_exp, solver_exp, res_exp, unb_exp, _ = _rollout_expert(
                x0,
                expert_u,
                target_seq,
                li,
                cfg_dyn,
            )
            if args.debug_first_step and debug_payload is not None:
                print(f"debug_first_step ({name}): {debug_payload}")

        pol_len = tip_pol.shape[0]
        exp_len = tip_exp.shape[0]
        target_pol = target_seq[:pol_len]
        target_exp = target_seq[:exp_len]

        pol_nonfinite = ~np.isfinite(tip_pol).all(axis=1) if pol_len else np.array([], dtype=bool)
        exp_nonfinite = ~np.isfinite(tip_exp).all(axis=1) if exp_len else np.array([], dtype=bool)
        pol_target_nonfinite = ~np.isfinite(target_pol).all(axis=1) if pol_len else np.array([], dtype=bool)
        exp_target_nonfinite = ~np.isfinite(target_exp).all(axis=1) if exp_len else np.array([], dtype=bool)

        pol_fail = (
            unb_pol.astype(bool)
            | pol_nonfinite
            | pol_target_nonfinite
            | (solver_pol != 0)
            | (res_pol > cfg_dyn.residual_threshold)
        )
        exp_fail = (
            unb_exp.astype(bool)
            | exp_nonfinite
            | exp_target_nonfinite
            | (solver_exp != 0)
            | (res_exp > cfg_dyn.residual_threshold)
        )

        pol_keep = ~pol_fail
        exp_keep = ~exp_fail

        err_pol = np.linalg.norm(tip_pol[pol_keep] - target_pol[pol_keep], axis=1) if pol_keep.any() else np.zeros((0,))
        err_exp = np.linalg.norm(tip_exp[exp_keep] - target_exp[exp_keep], axis=1) if exp_keep.any() else np.zeros((0,))

        mean_err = float(np.mean(err_pol)) if err_pol.size else float("nan")
        max_err = float(np.max(err_pol)) if err_pol.size else float("nan")
        base_mean = float(np.mean(err_exp)) if err_exp.size else float("nan")
        base_max = float(np.max(err_exp)) if err_exp.size else float("nan")
        policy_unbounded = float(np.mean(unb_pol)) if unb_pol.size else 0.0
        expert_unbounded = float(np.mean(unb_exp)) if unb_exp.size else 0.0

        if args.output_prefix:
            out_prefix = f"{args.output_prefix}_{name}"
        else:
            out_prefix = f"bc_eval_{name}_{stamp}"
        out_path = os.path.join("output_data", f"{out_prefix}.npz")
        np.savez(
            out_path,
            x_policy=x_pol,
            tip_policy=tip_pol,
            u_policy=u_pol,
            u_pred_raw=u_raw_pol,
            u_pred_used=u_pol,
            target_xyz=target_pol,
            solver_exit_policy=solver_pol,
            residual_norm_policy=res_pol,
            unbounded_detected_policy=unb_pol,
            policy_keep_mask=pol_keep,
            x_expert=x_exp,
            tip_expert=tip_exp,
            u_expert=expert_u[: tip_exp.shape[0]],
            solver_exit_expert=solver_exp,
            residual_norm_expert=res_exp,
            unbounded_detected_expert=unb_exp,
            expert_keep_mask=exp_keep,
            error_policy=err_pol,
            error_expert=err_exp,
            attempted_steps=int(target_seq.shape[0]),
            kept_steps_policy=int(pol_keep.sum()),
            dropped_steps_policy=int(pol_fail.sum()),
            kept_steps_expert=int(exp_keep.sum()),
            dropped_steps_expert=int(exp_fail.sum()),
            prev_u_init=prev_u_init,
            prev_u_source=prev_u_source,
            debug_first_step=debug_payload if debug_payload is not None else {},
            dt=dt,
            Li=li_val,
            start_idx=start_idx,
            subset_len=args.subset_len,
            ramp_len=ramp_len,
            umax=umax_safe,
            d_umax=d_umax_safe,
            dataset_name=os.path.basename(path),
            init_mode=args.init_mode,
            prefix_failed=prefix_failed if prefix_failed is not None else "",
            preroll_steps=preroll_steps,
            preroll_start_idx=preroll_start_idx,
            preroll_end_idx=preroll_end_idx,
            u_preroll=preroll_u,
        )

        if args.plot:
            plot_prefix = os.path.join("output_data", f"{out_prefix}")
            _plot("output_data", out_prefix, target_pol, tip_pol, tip_exp, err_pol, err_exp)
        else:
            plot_prefix = None

        if args.report_prefix:
            report_prefix = f"{args.report_prefix}_{name}"
        else:
            report_prefix = f"eval_bc_policy_{name}_{stamp}"
        report_path = os.path.join("docs", "control", "run_reports", f"{report_prefix}.md")
        write_run_report(
            report_path=report_path,
            command=" ".join(sys.argv),
            dataset_name=os.path.basename(path),
            start_idx=start_idx,
            end_idx=min(start_idx + args.subset_len, targets.shape[0]),
            n_total=targets.shape[0],
            li_mm=li_val,
            dt=dt,
            umax=umax_safe,
            d_umax=d_umax_safe,
            mean_error=mean_err,
            max_error=max_err,
            artifact_npz=out_path,
            plot_prefix=plot_prefix,
            init_from_ramp=args.init_mode == "init_from_ramp",
            ramp_len=ramp_len,
            max_wall_hit=False,
            baseline_mean=base_mean,
            baseline_max=base_max,
            extra_entries={
                "policy_unbounded_frac": policy_unbounded,
                "expert_unbounded_frac": expert_unbounded,
                "policy_steps_total": int(pol_len),
                "policy_steps_kept": int(pol_keep.sum()) if pol_len else 0,
                "policy_steps_dropped": int(pol_fail.sum()) if pol_len else 0,
                "expert_steps_total": int(exp_len),
                "expert_steps_kept": int(exp_keep.sum()) if exp_len else 0,
                "expert_steps_dropped": int(exp_fail.sum()) if exp_len else 0,
                "attempted_steps": int(target_seq.shape[0]),
                "prev_u_source": prev_u_source,
                "init_mode": args.init_mode,
                "base_start_idx": base_start_idx,
                "preroll_steps": preroll_steps,
                "preroll_start_idx": preroll_start_idx,
                "preroll_end_idx": preroll_end_idx,
                "debug_first_step": debug_payload if debug_payload is not None else "n/a",
                "eval_go": bool(pol_keep.sum() > 0),
                "model_path": args.model,
                "prefix_failed": prefix_failed if prefix_failed is not None else "",
                **(bc_counts or {}),
            },
        )

        print(f"{name}: saved {out_path}")
        print(f"{name}: report {report_path}")


if __name__ == "__main__":
    main()
