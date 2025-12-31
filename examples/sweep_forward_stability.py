import argparse
import csv
import io
import os
import sys
import time
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime

import numpy as np
import torch

from crm_diffsims.control.step_profiler import run_rollout_profile
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext


def _parse_list(text):
    return [float(item) for item in text.split(",") if item.strip()]


def _rate_limit(u_seq, d_umax):
    u_limited = u_seq.copy()
    for t in range(1, u_limited.shape[0]):
        delta = u_limited[t] - u_limited[t - 1]
        delta = np.clip(delta, -d_umax, d_umax)
        u_limited[t] = u_limited[t - 1] + delta
    return u_limited


def _build_u_seq(T, umax, d_umax, seed, dataset_currents=None):
    rng = np.random.default_rng(seed)
    if dataset_currents is not None:
        base = dataset_currents[:T].copy()
        if base.shape[0] < T:
            pad = np.repeat(base[-1:], T - base.shape[0], axis=0)
            base = np.concatenate([base, pad], axis=0)
        scale = umax / max(np.max(np.abs(base)), 1e-6)
        u = base * scale
    else:
        t = np.linspace(0.0, 2.0 * np.pi, T)
        u = np.zeros((T, 3), dtype=np.float64)
        u[:, 0] = umax * np.sin(t)
        u[:, 1] = 0.5 * umax * np.sin(2.0 * t + 0.5)
        u[:, 2] = 0.3 * umax * np.cos(1.5 * t)
        u += 0.01 * umax * rng.standard_normal(u.shape)
        u = np.clip(u, -umax, umax)

    u = _rate_limit(u, d_umax)
    return u


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, default=None)
    parser.add_argument("--umax_list", type=str, required=True)
    parser.add_argument("--d_umax_list", type=str, required=True)
    parser.add_argument("--T", type=int, default=50)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--single", action="store_true", help="Run a single umax/d_umax pair")
    parser.add_argument("--umax", type=float, default=None)
    parser.add_argument("--d_umax", type=float, default=None)
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, li = build_state()

    dataset_currents = None
    if args.dataset:
        data = np.load(args.dataset)
        if "currents" in data:
            dataset_currents = data["currents"]

    umax_list = _parse_list(args.umax_list)
    d_umax_list = _parse_list(args.d_umax_list)
    if args.single:
        if args.umax is None or args.d_umax is None:
            raise SystemExit("--single requires --umax and --d_umax")
        umax_list = [float(args.umax)]
        d_umax_list = [float(args.d_umax)]

    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    failure_dir = os.path.join(out_dir, "forward_failures")
    os.makedirs(failure_dir, exist_ok=True)

    csv_path = os.path.join(out_dir, "forward_stability_sweep.csv")
    fieldnames = [
        "umax",
        "d_umax",
        "T",
        "mean_step_time",
        "max_step_time",
        "nan_detected",
        "unbounded_detected",
        "fail_index",
    ]

    rows = []
    for umax in umax_list:
        for d_umax in d_umax_list:
            u_np = _build_u_seq(args.T, umax, d_umax, args.seed, dataset_currents)
            u_seq = torch.tensor(u_np, dtype=torch.float64).view(args.T, 1, 3)

            stdout_buf = io.StringIO()
            stderr_buf = io.StringIO()
            start = time.perf_counter()
            with redirect_stdout(stdout_buf), redirect_stderr(stderr_buf):
                result = run_rollout_profile(
                    x0,
                    u_seq,
                    li,
                    cfg_dyn,
                    label=f"umax{umax}_d{d_umax}",
                    save_dir=failure_dir,
                )
            elapsed = time.perf_counter() - start

            combined = stdout_buf.getvalue() + stderr_buf.getvalue()
            unbounded = "Unbounded" in combined

            row = {
                "umax": umax,
                "d_umax": d_umax,
                "T": args.T,
                "mean_step_time": result.mean_step_time,
                "max_step_time": result.max_step_time,
                "nan_detected": bool(result.nan_steps),
                "unbounded_detected": unbounded,
                "fail_index": -1 if result.fail_index is None else result.fail_index,
            }
            rows.append(row)

            if args.single:
                print(
                    f"single umax={umax:.3f} d_umax={d_umax:.3f} "
                    f"mean_step={result.mean_step_time:.4f}s max_step={result.max_step_time:.4f}s "
                    f"nan={bool(result.nan_steps)} unbounded={unbounded} elapsed={elapsed:.2f}s"
                )
            else:
                print(
                    f"Sweep umax={umax:.3f} d_umax={d_umax:.3f} "
                    f"mean_step={result.mean_step_time:.4f}s max_step={result.max_step_time:.4f}s "
                    f"nan={bool(result.nan_steps)} unbounded={unbounded} elapsed={elapsed:.2f}s"
                )

            if result.fail_index is None and not unbounded:
                json_path = os.path.join(failure_dir, f"flight_umax{umax}_d{d_umax}.json")
                npz_path = os.path.join(failure_dir, f"flight_umax{umax}_d{d_umax}.npz")
                if os.path.exists(json_path):
                    os.remove(json_path)
                if os.path.exists(npz_path):
                    os.remove(npz_path)

    if not args.single:
        with open(csv_path, "w", newline="", encoding="ascii") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if not args.single:
        print(f"Saved sweep results to {csv_path} ({stamp})")

    if args.single and any(r["nan_detected"] for r in rows):
        sys.exit(1)


if __name__ == "__main__":
    main()
