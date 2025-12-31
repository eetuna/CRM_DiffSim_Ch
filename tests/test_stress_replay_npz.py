import glob
import csv
from pathlib import Path

import numpy as np
import torch
from torch.utils.cpp_extension import load

from crm_diffsims.dynamics.step import build_stable_state, crm_step, load_dyn_ext


def _parse_param_file(path):
    values = {}
    with open(path, "r", encoding="ascii") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            key = parts[0]
            vals = []
            for item in parts[1:]:
                try:
                    vals.append(float(item))
                except ValueError:
                    vals = []
                    break
            if vals:
                values[key] = vals
    return values


def _act_inertia_from_params(params):
    outer_r = params["oRlist"][0]
    inner_r = params["iRlist"][0]
    act_mass = params["ActMass"][0]
    seg_len = params["SegmentLengths"][1]
    i_zz = 0.5 * act_mass * (outer_r * outer_r + inner_r * inner_r)
    i_xx = 0.25 * act_mass * (outer_r * outer_r + inner_r * inner_r) + (1.0 / 12.0) * act_mass * seg_len * seg_len
    return [
        i_xx, 0.0, 0.0,
        0.0, i_xx, 0.0,
        0.0, 0.0, i_zz,
    ]


def _build_cfg(ext, integration_step, delta_t):
    repo_root = Path(__file__).resolve().parents[1]
    cath_params_path = repo_root / "catheterdata" / "CatheterParameterSet_1_dyn.txt"
    cath_config_path = repo_root / "catheterdata" / "CatheterSpatialConfiguration_1.txt"

    params = _parse_param_file(cath_params_path)
    act_inertia = _act_inertia_from_params(params)

    damping = [
        12.1761626666366, 12.1761626666366, 284.429938756989,
        0.0304776127617393, 0.0304776127617393, 0.00502712804532508,
    ]

    return ext.DynamicsConfig(
        str(cath_params_path),
        str(cath_config_path),
        float(integration_step),
        [0.0, 0.0, 0.0],
        act_inertia,
        damping,
        float(delta_t),
        1e-2,
        1e-5,
    )


def _build_equilibrium_cfg(ext, integration_step):
    repo_root = Path(__file__).resolve().parents[1]
    cath_params_path = repo_root / "catheterdata" / "CatheterParameterSet_1_dyn.txt"
    cath_config_path = repo_root / "catheterdata" / "CatheterSpatialConfiguration_1.txt"
    return ext.EquilibriumConfig(
        str(cath_params_path),
        str(cath_config_path),
        float(integration_step),
        [0.0, 0.0, 0.0],
        1e-5,
    )


def _equilibrium_seed_state(ext, cfg_eq, li, u0):
    u_in = torch.tensor(u0, dtype=torch.float64).view(1, 1, 3)
    p_tip, deltau0, ftip, localmin, residual_norm = ext.crm_equilibrium_forward(u_in, li, cfg_eq)
    x0 = torch.zeros((1, 24 * 1 + 15), dtype=torch.float64)
    x0[0, 12:15] = p_tip[0]
    x0[0, 15:24] = torch.tensor([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=torch.float64)
    x0[0, 24:27] = p_tip[0]
    x0[0, 27:36] = torch.tensor([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=torch.float64)
    return x0


def _load_equilibrium_ext():
    repo_root = Path(__file__).resolve().parents[1]
    build_dir = repo_root / "build_torch_ext_equil"
    build_dir.mkdir(parents=True, exist_ok=True)
    sources = [
        str(repo_root / "src" / "CRM_TorchEquilibrium.cpp"),
        str(repo_root / "src" / "CRM_IVPSolver.cpp"),
        str(repo_root / "src" / "CRM_BVPSolver.cpp"),
        str(repo_root / "src" / "CRM_IVPJacobian.cpp"),
        str(repo_root / "src" / "CRM_ForwardKinematics.cpp"),
        str(repo_root / "src" / "CRM_CatheterClass.cpp"),
        str(repo_root / "src" / "CRM_SupportFunctions.cpp"),
        str(repo_root / "src" / "CoilDynamics_Defs.cpp"),
        str(repo_root / "numerical" / "minpack.cpp"),
        str(repo_root / "numerical" / "minpack_DYN_Defs.cpp"),
    ]
    include_dirs = [
        str(repo_root / "src"),
        str(repo_root / "numerical"),
    ]
    return load(
        name="crm_equilibrium_ext",
        sources=sources,
        extra_include_paths=include_dirs,
        extra_cflags=["-O3", "-std=c++17"],
        build_directory=str(build_dir),
        verbose=False,
    )

def _load_workspace_nn(repo_root):
    files = sorted((repo_root / "data").glob("workspace_fk_ins94.3_*.npz"))
    if not files:
        return None
    data = np.load(files[0])
    return {
        "U": data["U"],
        "P": data["P"],
        "conv": data["conv"],
        "name": files[0].name,
    }


def _nearest_neighbor(u_vec, U):
    diff = U - u_vec[None, :]
    dists = np.sum(diff * diff, axis=1)
    return int(np.argmin(dists))


def _apply_workspace_hint(x_base, u_vec, workspace):
    x_in = x_base.clone()
    if workspace is None:
        return x_in
    idx = _nearest_neighbor(u_vec, workspace["U"])
    p = workspace["P"][idx]
    x_in[0, 24:27] = torch.tensor(p, dtype=x_in.dtype)
    return x_in


def _rollout(ext, cfg, li, x_base, u_seq, tip_dyn, workspace, init_mode, use_crm_step):
    residuals = []
    exit_flags = []
    errors = []
    tip_preds = []

    x_t = x_base
    for t in range(u_seq.shape[0]):
        u_t = u_seq[t:t + 1]
        if init_mode == "zero":
            x_in = x_base.clone() if t == 0 else x_t
        elif init_mode == "workspace":
            x_in = _apply_workspace_hint(x_base, u_seq[t].numpy(), workspace) if t == 0 else x_t
        elif init_mode == "equilibrium_seed":
            x_in = x_base.clone() if t == 0 else x_t
        else:
            x_in = x_t

        x_in = x_in.clone().detach()
        if init_mode != "warm":
            x_in[0, 24:27] = torch.tensor(tip_dyn[t], dtype=x_in.dtype)
        if use_crm_step:
            x_tp1 = crm_step(x_in, u_t, li, cfg)
            x_tp1_ref, eta, exit_flag, residual_norm, u0, tau = ext.crm_step_forward(x_in, u_t, li, cfg)
        else:
            x_tp1, eta, exit_flag, residual_norm, u0, tau = ext.crm_step_forward(x_in, u_t, li, cfg)
        exit_flags.append(int(exit_flag.item()))
        residuals.append(float(residual_norm.item()))

        tip_pred = x_tp1.detach()[0, 24:27].cpu().numpy()
        tip_ref = tip_dyn[t]
        errors.append(float(np.linalg.norm(tip_pred - tip_ref)))
        tip_preds.append(tip_pred.copy())

        x_t = x_tp1.detach()

    return {
        "exit_flags": np.array(exit_flags),
        "residuals": np.array(residuals),
        "errors": np.array(errors),
        "tip_preds": np.array(tip_preds),
    }


def test_stress_replay_npz():
    repo_root = Path(__file__).resolve().parents[1]
    data_files = sorted(glob.glob(str(repo_root / "data" / "dyn_fk_*.npz")))
    assert data_files, "No dyn_fk_*.npz files found in ./data"

    ext = load_dyn_ext()
    ext_eq = _load_equilibrium_ext()
    workspace = _load_workspace_nn(repo_root)
    out_dir = repo_root / "output_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "stress_outliers_943.csv"
    out_rows = []

    for fname in data_files:
        data = np.load(fname)
        u_seq = torch.tensor(data["currents"], dtype=torch.float64).view(-1, 1, 3)
        tip_dyn = data["tip_dyn"]
        t_arr = data["t"] if "t" in data.files else np.arange(u_seq.shape[0], dtype=np.float64)
        assert np.isfinite(tip_dyn).all(), f"tip_dyn has NaN/Inf in {Path(fname).name}"

        li = torch.tensor([float(np.array(data["insertion_length"]).item())], dtype=torch.float64)
        dt = float(np.array(data["dt"]).item()) if "dt" in data.files else 0.05
        integration_step = float(np.array(data["integration_step_size"]).item()) if "integration_step_size" in data.files else 0.2
        cfg = _build_cfg(ext, integration_step, dt)
        cfg_eq = _build_equilibrium_cfg(ext_eq, integration_step)

        x_base, _ = build_stable_state(ext, cfg)
        x_base_hint = x_base.clone()
        x_base_hint[0, 24:27] = torch.tensor(tip_dyn[0], dtype=x_base_hint.dtype)
        x_eq = _equilibrium_seed_state(ext_eq, cfg_eq, li, u_seq[0, 0].cpu().numpy())

        stats = {}
        for mode in ("zero", "warm", "workspace", "equilibrium_seed"):
            if mode == "equilibrium_seed":
                base = x_eq
            elif mode == "warm":
                base = x_base
            else:
                base = x_base_hint
            use_crm = mode in ("warm", "equilibrium_seed")
            stats[mode] = _rollout(ext, cfg, li, base, u_seq, tip_dyn, workspace, mode, use_crm)

        def summarize(label, s):
            exit_flags = s["exit_flags"]
            residuals = s["residuals"]
            errors = s["errors"]
            print(
                f"{Path(fname).name} {label}: "
                f"conv_rate={np.mean(exit_flags == 0):.3f} "
                f"max_res={np.max(residuals):.3e} "
                f"mean_err={np.mean(errors):.3e} "
                f"max_err={np.max(errors):.3e}"
            )

        for mode, s in stats.items():
            summarize(mode, s)

        best_mode = max(
            stats.keys(),
            key=lambda m: (np.mean(stats[m]["exit_flags"] == 0), -np.mean(stats[m]["residuals"])),
        )
        print(f"{Path(fname).name} best_init={best_mode}")

        for mode, s in stats.items():
            bad = (s["exit_flags"] != 0) & (s["residuals"] > cfg.residual_threshold * 5.0)
            assert not np.any(bad), f"non-convergence in {Path(fname).name} mode={mode}"

        residual_threshold = max(cfg.residual_threshold * 5.0, 5e-2)
        for mode, s in stats.items():
            assert np.max(s["residuals"]) <= residual_threshold, (
                f"residual spike in {Path(fname).name} mode={mode}"
            )

        catastrophic = 50.0
        for mode, s in stats.items():
            assert np.max(s["errors"]) <= catastrophic, (
                f"tip divergence in {Path(fname).name} mode={mode}"
            )

        k = 20
        eq_errors = stats["equilibrium_seed"]["errors"]
        for mode, s in stats.items():
            idxs = np.argsort(s["errors"])[-k:][::-1]
            for i in idxs:
                tip_err = float(s["errors"][i])
                res_norm = float(s["residuals"][i])
                exit_flag = int(s["exit_flags"][i])
                nonconv = exit_flag != 0 and res_norm > cfg.residual_threshold * 5.0
                inconsistency = (res_norm <= cfg.residual_threshold * 5.0) and (tip_err > 10.0)
                eq_drop = (i <= 5) and (eq_errors[i] <= max(1.0, 0.2 * tip_err))
                if inconsistency and eq_drop and mode != "equilibrium_seed":
                    classification = "expected_initial_state_mismatch"
                elif nonconv:
                    classification = "expected_nonconvergence"
                elif inconsistency:
                    classification = "inconsistency"
                else:
                    classification = "ok"

                u_prev = u_seq[i - 1, 0].cpu().numpy().tolist() if i > 0 else [float("nan")] * 3
                u_curr = u_seq[i, 0].cpu().numpy().tolist()
                u_next = u_seq[i + 1, 0].cpu().numpy().tolist() if i + 1 < u_seq.shape[0] else [float("nan")] * 3
                out_rows.append({
                    "dataset_name": Path(fname).name,
                    "init_mode": mode,
                    "index": int(i),
                    "time": float(t_arr[i]),
                    "tip_error": tip_err,
                    "tip_dyn_saved": tip_dyn[i].tolist(),
                    "tip_dyn_replay": s["tip_preds"][i].tolist(),
                    "solver_exit": exit_flag,
                    "residual_norm": res_norm,
                    "u_prev": u_prev,
                    "u_curr": u_curr,
                    "u_next": u_next,
                    "classification": classification,
                })

    inconsistencies = [r for r in out_rows if r["classification"] == "inconsistency"]
    mismatch = [r for r in out_rows if r["classification"] == "expected_initial_state_mismatch"]
    nonconv = [r for r in out_rows if r["classification"] == "expected_nonconvergence"]
    by_mode = {}
    for r in inconsistencies:
        by_mode[r["init_mode"]] = by_mode.get(r["init_mode"], 0) + 1
    print(f"inconsistencies_detected={len(inconsistencies)}")
    print(f"mismatch_detected={len(mismatch)}")
    print(f"nonconvergence_detected={len(nonconv)}")
    print(f"inconsistencies_by_mode={by_mode}")

    with open(out_path, "w", encoding="ascii", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dataset_name",
                "init_mode",
                "index",
                "time",
                "tip_error",
                "tip_dyn_saved",
                "tip_dyn_replay",
                "solver_exit",
                "residual_norm",
                "u_prev",
                "u_curr",
                "u_next",
                "classification",
            ],
        )
        writer.writeheader()
        writer.writerows(out_rows)
