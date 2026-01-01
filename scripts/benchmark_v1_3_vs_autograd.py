import json
import time
from pathlib import Path

import numpy as np
import torch

from crm_diffsims.dynamics.linearize_v1_3 import linearize_v1_3
from crm_diffsims.dynamics.step import build_cfg, build_stable_state, load_dyn_ext
from crm_diffsims.dynamics.step import crm_step


def _time_block(fn, n):
    times = []
    for _ in range(n):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return times


def main():
    torch.manual_seed(0)
    np.random.seed(0)

    ext = load_dyn_ext()
    cfg = build_cfg(ext)
    x0, li = build_stable_state(ext, cfg)

    u = torch.zeros((1, 1, 3), dtype=torch.float64)
    u[0, 0, 2] = 0.005

    def autograd_jac():
        x_t = x0.detach().requires_grad_(True)
        u_t = u.detach().requires_grad_(True)

        def _dyn_x(x_in):
            return crm_step(x_in, u_t, li, cfg)

        def _dyn_u(u_in):
            return crm_step(x_t, u_in, li, cfg)

        try:
            a = torch.autograd.functional.jacobian(_dyn_x, x_t, create_graph=False, vectorize=True)
            b = torch.autograd.functional.jacobian(_dyn_u, u_t, create_graph=False, vectorize=True)
        except RuntimeError:
            a = torch.autograd.functional.jacobian(_dyn_x, x_t, create_graph=False, vectorize=False)
            b = torch.autograd.functional.jacobian(_dyn_u, u_t, create_graph=False, vectorize=False)
        return a, b

    def v1_3_op_apply():
        a_op, b_op = linearize_v1_3(x0, u, li, cfg, return_dense=False)
        dx = torch.randn_like(x0) * 1e-3
        du = torch.randn_like(u) * 1e-3
        _ = a_op.apply(dx)
        _ = b_op.apply(du)

    # Warmup
    autograd_jac()
    v1_3_op_apply()

    n = 1
    auto_times = _time_block(autograd_jac, n)
    v1_times = _time_block(v1_3_op_apply, n)

    auto_mean = float(np.mean(auto_times))
    v1_mean = float(np.mean(v1_times))
    speedup = auto_mean / v1_mean if v1_mean > 0 else float("inf")

    result = {
        "autograd": {"mean_sec": auto_mean, "min_sec": float(np.min(auto_times)), "max_sec": float(np.max(auto_times))},
        "v1_3_operator": {
            "mean_sec": v1_mean,
            "min_sec": float(np.min(v1_times)),
            "max_sec": float(np.max(v1_times)),
        },
        "speedup": speedup,
        "iterations": n,
    }

    out_dir = Path("output_data")
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "v1_3_benchmark.json"
    txt_path = out_dir / "v1_3_benchmark.txt"

    json_path.write_text(json.dumps(result, indent=2), encoding="ascii")
    txt_path.write_text(
        "v1.3 vs autograd benchmark\n"
        f"autograd mean: {auto_mean:.6f}s (min {min(auto_times):.6f}s, max {max(auto_times):.6f}s)\n"
        f"v1.3 operator mean: {v1_mean:.6f}s (min {min(v1_times):.6f}s, max {max(v1_times):.6f}s)\n"
        f"speedup: {speedup:.1f}x\n",
        encoding="ascii",
    )

    print(txt_path.read_text(encoding="ascii"))


if __name__ == "__main__":
    main()
