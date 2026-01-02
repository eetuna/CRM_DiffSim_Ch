import argparse
import glob
import os
import sys
from datetime import datetime

if __package__ is None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    python_root = os.path.join(repo_root, "python")
    if python_root not in sys.path:
        sys.path.insert(0, python_root)

import numpy as np

from crm_diffsims.control.run_report import write_run_report


def _infer_dataset_label(path, payload):
    name = str(payload.get("dataset_name", ""))
    if name:
        return os.path.basename(name)
    base = os.path.basename(path)
    if "circle" in base:
        return "circle"
    if "lemniscate" in base:
        return "lemniscate"
    return base


def _drop_reasons(payload, residual_threshold=0.01):
    keep = payload.get("policy_keep_mask")
    if keep is None:
        return {}, 0, 0, 0
    keep = keep.astype(bool)
    total = keep.shape[0]
    kept = int(keep.sum())
    dropped = int(total - kept)
    if dropped == 0:
        return {
            "unbounded": 0,
            "solver_exit": 0,
            "residual": 0,
            "nonfinite_tip": 0,
            "target_nonfinite": 0,
        }, total, kept, dropped

    drop_mask = ~keep
    unbounded = payload.get("unbounded_detected_policy", np.zeros((total,), dtype=bool)).astype(bool)
    solver = payload.get("solver_exit_policy", np.zeros((total,), dtype=np.int32))
    residual = payload.get("residual_norm_policy", np.zeros((total,), dtype=np.float64))
    tip = payload.get("tip_policy", np.zeros((total, 3), dtype=np.float64))
    target = payload.get("target_xyz", np.zeros((total, 3), dtype=np.float64))

    reasons = {
        "unbounded": int((unbounded & drop_mask).sum()),
        "solver_exit": int(((solver != 0) & drop_mask).sum()),
        "residual": int(((residual > residual_threshold) & drop_mask).sum()),
        "nonfinite_tip": int((~np.isfinite(tip).all(axis=1) & drop_mask).sum()),
        "target_nonfinite": int((~np.isfinite(target).all(axis=1) & drop_mask).sum()),
    }
    return reasons, total, kept, dropped


def _format_row(row):
    return (
        f"{row['dataset']:<12} "
        f"{row['start_idx']:<8} {row['end_idx']:<8} "
        f"{row['attempted_steps']:<8} {row['kept_steps']:<8} {row['kept_frac']:<8} "
        f"{row['mean_error_kept']:<14} {row['max_error_kept']:<14} "
        f"{row['drop_reasons_counts']} "
        f"{row['policy_unbounded_frac']:<8} {row['expert_unbounded_frac']:<8}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+", help="Eval NPZ paths (glob allowed)")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    paths = []
    for pattern in args.paths:
        matches = glob.glob(pattern)
        if matches:
            paths.extend(matches)
        else:
            paths.append(pattern)

    rows = []
    for path in sorted(set(paths)):
        if not os.path.exists(path):
            print(f"missing: {path}")
            continue
        payload = np.load(path, allow_pickle=True)
        dataset = _infer_dataset_label(path, payload)
        start_idx = int(payload.get("start_idx", 0))
        subset_len = int(payload.get("subset_len", 0))
        end_idx = start_idx + subset_len
        attempted_steps = int(payload.get("attempted_steps", payload.get("target_xyz", np.zeros((0,))).shape[0]))
        kept_steps = int(payload.get("kept_steps_policy", payload.get("policy_keep_mask", np.zeros((0,))).sum()))
        kept_frac = kept_steps / max(attempted_steps, 1)
        mean_err = payload.get("error_policy", np.array([]))
        max_err = mean_err
        mean_error_kept = float(np.mean(mean_err)) if mean_err.size else float("nan")
        max_error_kept = float(np.max(max_err)) if max_err.size else float("nan")
        reasons, _, _, _ = _drop_reasons(payload)
        policy_unb = payload.get("unbounded_detected_policy", np.array([]))
        expert_unb = payload.get("unbounded_detected_expert", np.array([]))
        policy_unb_frac = float(np.mean(policy_unb)) if policy_unb.size else 0.0
        expert_unb_frac = float(np.mean(expert_unb)) if expert_unb.size else 0.0

        rows.append(
            {
                "dataset": dataset,
                "start_idx": start_idx,
                "end_idx": end_idx,
                "attempted_steps": attempted_steps,
                "kept_steps": kept_steps,
                "kept_frac": f"{kept_frac:.3f}",
                "mean_error_kept": f"{mean_error_kept:.3f}",
                "max_error_kept": f"{max_error_kept:.3f}",
                "drop_reasons_counts": reasons,
                "policy_unbounded_frac": f"{policy_unb_frac:.3f}",
                "expert_unbounded_frac": f"{expert_unb_frac:.3f}",
            }
        )

    header = (
        "dataset      start    end      attempted kept     kept_frac "
        "mean_error_kept max_error_kept drop_reasons_counts policy_unb expert_unb"
    )
    print(header)
    for row in rows:
        print(_format_row(row))

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join("docs", "control", "run_reports", f"summary_{stamp}.md")
    os.makedirs(os.path.dirname(report_path), exist_ok=True)

    npz_path = os.path.join("output_data", f"summary_{stamp}.npz")
    np.savez(npz_path, rows=np.array(rows, dtype=object))

    write_run_report(
        report_path=report_path,
        command=" ".join(sys.argv),
        dataset_name="bc_eval_summary",
        start_idx=0,
        end_idx=0,
        n_total=len(rows),
        li_mm=0.0,
        dt=0.0,
        umax=0.0,
        d_umax=0.0,
        mean_error=0.0,
        max_error=0.0,
        artifact_npz=npz_path,
        plot_prefix=None,
        init_from_ramp=False,
        ramp_len=0,
        max_wall_hit=False,
        extra_entries={
            "summary_rows": len(rows),
        },
    )

    lines = ["# BC Eval Summary", "", header]
    lines.extend(_format_row(row) for row in rows)
    lines.append("")
    with open(report_path, "a", encoding="ascii") as f:
        f.write("\n".join(lines))
    print(f"Summary report: {report_path}")

    if args.plot and rows:
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available; skipping plots")
            return
        start_vals = [r["start_idx"] for r in rows]
        kept_vals = [float(r["kept_frac"]) for r in rows]
        mean_vals = [float(r["mean_error_kept"]) for r in rows]

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(start_vals, kept_vals, marker="o")
        ax.set_title("Kept fraction vs start_idx")
        ax.set_xlabel("start_idx")
        ax.set_ylabel("kept_frac")
        fig.tight_layout()
        out_prefix = os.path.join("output_data", f"summary_kept_frac_{stamp}")
        fig.savefig(out_prefix + ".png", dpi=300)
        fig.savefig(out_prefix + ".pdf")

        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(start_vals, mean_vals, marker="o")
        ax.set_title("Mean error vs start_idx")
        ax.set_xlabel("start_idx")
        ax.set_ylabel("mean_error_kept")
        fig.tight_layout()
        out_prefix = os.path.join("output_data", f"summary_mean_error_{stamp}")
        fig.savefig(out_prefix + ".png", dpi=300)
        fig.savefig(out_prefix + ".pdf")


if __name__ == "__main__":
    main()
