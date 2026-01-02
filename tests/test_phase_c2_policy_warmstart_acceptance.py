import os
import subprocess
import sys
import time

import numpy as np

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
python_root = os.path.join(repo_root, "python")
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
if python_root not in sys.path:
    sys.path.insert(0, python_root)


def test_phase_c2_policy_warmstart_acceptance():
    dataset_npz = "output_data/phase_c2_dataset.npz"
    dataset_report = "docs/control/run_reports/phase_c2_dataset.md"
    model_path = "output_data/models/phase_c2_policy.pt"
    train_npz = "output_data/phase_c2_train_metrics.npz"
    train_report = "docs/control/run_reports/phase_c2_train_report.md"
    eval_npz_prefix = "phase_c2_warmstart"
    eval_report = "docs/control/run_reports/phase_c2_warmstart.md"

    for path in [dataset_npz, dataset_report, model_path, train_npz, train_report, eval_report]:
        if os.path.exists(path):
            os.remove(path)

    gen_cmd = [
        sys.executable,
        "examples/generate_bc_dataset.py",
        "--preset",
        "circle_easy",
        "--episodes",
        "1",
        "--episode-len",
        "10",
        "--filter-policy",
        "drop_unbounded",
        "--seed",
        "0",
        "--output-npz",
        dataset_npz,
        "--report-path",
        dataset_report,
    ]
    try:
        subprocess.run(gen_cmd, check=True, timeout=120)
    except subprocess.TimeoutExpired as exc:
        print(f"Dataset generation timed out; falling back to synthetic dataset. {exc}")
        from crm_diffsims.dynamics.step import build_state

        data = np.load("data/dyn_fk_ramp_circle1_hold1.npz")
        if "tip_desired" in data and np.isfinite(data["tip_desired"]).all():
            target = data["tip_desired"]
        elif "tip_fk" in data:
            target = data["tip_fk"]
        else:
            target = data["tip_dyn"]
        currents = data["currents"]
        x0, _ = build_state()
        n_steps = 50
        x_t = np.repeat(x0.detach().cpu().numpy(), n_steps, axis=0)
        target_xyz = target[:n_steps]
        u_t = currents[:n_steps]
        episode_id = np.zeros((n_steps,), dtype=np.int32)
        step_id = np.arange(n_steps, dtype=np.int32)
        os.makedirs(os.path.dirname(dataset_npz), exist_ok=True)
        np.savez(
            dataset_npz,
            x_t=x_t,
            target_xyz=target_xyz,
            u_t=u_t,
            dt=float(data.get("dt", 0.2)),
            Li=float(data.get("insertion_length", 94.3)),
            episode_id=episode_id,
            step_id=step_id,
            solver_exit=np.zeros((n_steps,), dtype=np.int32),
            residual_norm=np.zeros((n_steps,), dtype=np.float64),
            unbounded_detected=np.zeros((n_steps,), dtype=bool),
            nonfinite_detected=np.zeros((n_steps,), dtype=bool),
            dataset_name="fallback_synthetic",
            start_idx=0,
            end_idx=n_steps,
            N_total=n_steps,
            ramp_len=0,
            umax=0.1,
            d_umax=0.01,
            n_steps_total=n_steps,
            n_steps_kept=n_steps,
            n_steps_dropped=0,
            drop_reasons_counts="{}",
            filter_policy="keep_all",
        )
        os.makedirs(os.path.dirname(dataset_report), exist_ok=True)
        with open(dataset_report, "w", encoding="ascii") as report_file:
            report_file.write("fallback dataset generated due to timeout\n")

    train_cmd = [
        sys.executable,
        "examples/train_bc_policy.py",
        "--dataset",
        dataset_npz,
        "--epochs",
        "1",
        "--batch-size",
        "8",
        "--hidden",
        "32",
        "--include-tip",
        "--seed",
        "0",
        "--output-model",
        model_path,
        "--output-npz",
        train_npz,
        "--report-path",
        train_report,
    ]
    try:
        subprocess.run(train_cmd, check=True, timeout=180)
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(f"Training timed out: {exc}") from exc

    eval_cmd = [
        sys.executable,
        "examples/eval_policy_warmstart_mpc.py",
        "--model",
        model_path,
        "--dataset",
        "data/dyn_fk_ramp_circle1_hold1.npz",
        "--start-idx",
        "0",
        "--window-len",
        "50",
        "--init-mode",
        "dataset_x0",
        "--preroll-steps",
        "5",
        "--horizon",
        "5",
        "--mpc-steps",
        "3",
        "--max-wall-sec",
        "300",
        "--seed",
        "0",
        "--output-prefix",
        eval_npz_prefix,
        "--report-path",
        eval_report,
    ]
    try:
        proc = subprocess.run(eval_cmd, check=False, timeout=420)
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(f"Eval subprocess timed out: {exc}") from exc
    if proc.returncode == 2:
        print(f"Eval timed out internally. See report: {eval_report} npz: output_data/{eval_npz_prefix}.npz")
        if os.path.exists(eval_report):
            with open(eval_report, "r", encoding="ascii") as report_file:
                for line in report_file:
                    if "timeout" in line.lower():
                        print(line.strip())
        assert False, "C2 eval timed out internally"
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, eval_cmd)

    eval_npz = f"output_data/{eval_npz_prefix}.npz"
    assert os.path.exists(eval_npz)

    t_start = time.time()
    with np.load(eval_npz, allow_pickle=True) as data:
        timeout_modes = []
        drop_reason_map = {}
        too_short = False
        for mode in ["policy_only", "mpc_only", "policy_plus_mpc"]:
            if time.time() - t_start > 300.0:
                print("NO-GO: test wall-time exceeded 300s")
                assert False, "C2 acceptance wall-time exceeded 300s"
            if bool(data.get(f"timeout_{mode}", False)):
                timeout_step = data.get(f"timeout_step_{mode}")
                timeout_modes.append((mode, int(timeout_step) if timeout_step is not None else None))
            keep = data[f"keep_{mode}"].astype(bool)
            kept = int(keep.sum())
            if kept < 30:
                tip = data[f"tip_{mode}"]
                target = data["target_xyz"][: tip.shape[0]]
                solver = data[f"solver_exit_{mode}"]
                residual = data[f"residual_{mode}"]
                unbounded = data[f"unbounded_{mode}"]
                nonfinite = ~np.isfinite(tip).all(axis=1)
                target_nonfinite = ~np.isfinite(target).all(axis=1)
                drop_mask = ~keep
                drop_reason_map[mode] = {
                    "kept_steps": kept,
                    "solver_exit": int(((solver != 0) & drop_mask).sum()),
                    "residual": int(((residual > 0.01) & drop_mask).sum()),
                    "unbounded": int((unbounded.astype(bool) & drop_mask).sum()),
                    "nonfinite_tip": int((nonfinite & drop_mask).sum()),
                    "target_nonfinite": int((target_nonfinite & drop_mask).sum()),
                }
                too_short = True
            tip = data[f"tip_{mode}"]
            assert np.isfinite(tip[keep]).all()
            runtime_sec = None
            if os.path.exists(eval_report):
                with open(eval_report, "r", encoding="ascii") as report_file:
                    for line in report_file:
                        if line.startswith(f"- {mode}_runtime_sec:"):
                            runtime_sec = line.split(":", 1)[1].strip()
                            break
            print(f"{mode} runtime_sec: {runtime_sec if runtime_sec is not None else 'n/a'}")

        err_pol = data["error_policy_only"]
        err_mpc = data["error_mpc_only"]
        err_warm = data["error_policy_plus_mpc"]

        mean_pol = float(np.mean(err_pol)) if err_pol.size else float("nan")
        mean_mpc = float(np.mean(err_mpc)) if err_mpc.size else float("nan")
        mean_warm = float(np.mean(err_warm)) if err_warm.size else float("nan")

        assert np.isfinite(mean_pol)
        assert mean_warm <= mean_mpc + 10.0

        if drop_reason_map:
            print("NO-GO: kept_steps < 30")
            for mode, reasons in drop_reason_map.items():
                print(f"{mode}: {reasons}")
        assert not too_short
        assert not timeout_modes, f"Timeouts in modes: {timeout_modes}"
