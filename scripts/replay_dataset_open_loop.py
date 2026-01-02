import argparse
import os
import sys
import time
from datetime import datetime

import numpy as np
import torch

if __package__ is None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    python_root = os.path.join(repo_root, "python")
    if python_root not in sys.path:
        sys.path.insert(0, python_root)

from crm_diffsims.control.ilqr import tip_position
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


def _atomic_save_npz(path, **payload):
    tmp_path = f"{path}.tmp.npz"
    np.savez(tmp_path, **payload)
    os.replace(tmp_path, path)


def _maybe_plot(out_dir, stem, tip_pred, tip_dyn, errors):
    try:
        from crm_diffsims.control import plot_utils
    except Exception:
        print("matplotlib not installed; skipping plots")
        return
    plot_utils.save_tip_trajectory_plot(
        tip_pred,
        tip_dyn,
        title="Tip trajectory (replay)",
        out_path=f"{out_dir}/replay_{stem}_tip",
    )
    plot_utils.save_error_plot(
        errors,
        title="Tip tracking error (replay)",
        out_path=f"{out_dir}/replay_{stem}_error",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--npz", required=True, help="Dataset npz with currents and tip_dyn")
    parser.add_argument("--plot", action="store_true", help="Save plots to output_data/")
    args = parser.parse_args()

    torch.manual_seed(0)
    np.random.seed(0)

    data = np.load(args.npz)
    if "currents" not in data:
        raise RuntimeError("Dataset missing currents")
    if "tip_dyn" not in data:
        raise RuntimeError("Dataset missing tip_dyn")

    currents = data["currents"]
    tip_dyn = data["tip_dyn"]
    dt = float(data.get("dt", 0.2))
    li_val = float(data.get("insertion_length", 94.3))

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, _ = build_state()
    x0 = x0.clone()
    x0[0, 24:27] = torch.tensor(tip_dyn[0], dtype=x0.dtype)
    li = torch.tensor([li_val], dtype=torch.float64)

    n_steps = currents.shape[0]
    u_seq = torch.tensor(currents, dtype=torch.float64).view(n_steps, 1, 3)

    x_roll = [x0]
    tip_preds = [tip_position(x0)[0].detach().cpu().numpy()]
    solver_exit = []
    residual_norm = []
    start = time.time()
    for t in range(n_steps):
        x_tp1, cache = crm_step_v1_3_forward(x_roll[-1], u_seq[t : t + 1], li, cfg_dyn)
        solver_exit.append(int(cache["solver_exit"].item()))
        residual_norm.append(float(cache["residual_norm"].item()))
        x_roll.append(x_tp1)
        tip_preds.append(tip_position(x_tp1)[0].detach().cpu().numpy())
        if solver_exit[-1] != 0 or residual_norm[-1] > cfg_dyn.residual_threshold:
            print(f"Nonconvergence at step {t}: exit={solver_exit[-1]} residual={residual_norm[-1]:.3e}")
            break

    tip_preds = np.asarray(tip_preds)
    tip_dyn_used = tip_dyn[: tip_preds.shape[0]]
    errors = np.linalg.norm(tip_preds - tip_dyn_used, axis=1)
    print(f"Replay mean_err={errors.mean():.3e} max_err={errors.max():.3e} in {time.time() - start:.1f}s")

    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = f"{os.path.splitext(os.path.basename(args.npz))[0]}_{stamp}"
    out_path = f"{out_dir}/replay_{stem}.npz"
    _atomic_save_npz(
        out_path,
        tip_dyn=tip_dyn_used,
        tip_pred=tip_preds,
        errors=errors,
        currents=currents[: tip_preds.shape[0] - 1],
        dt=dt,
        Li_mm=li_val,
        solver_exit=np.array(solver_exit, dtype=np.int32),
        residual_norm=np.array(residual_norm, dtype=np.float64),
    )

    if args.plot:
        _maybe_plot(out_dir, stem, tip_preds, tip_dyn_used, errors)
        print(f"Saved plots to {out_dir}/replay_{stem}_*.{{png,pdf}}")


if __name__ == "__main__":
    main()
