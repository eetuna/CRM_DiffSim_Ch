import csv
import glob
from pathlib import Path

import numpy as np
import torch
from torch.utils.cpp_extension import load

from crm_diffsims.dynamics.step import build_state, crm_step, load_dyn_ext


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


def _sign(x, eps):
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def _workspace_sign_agreement(repo_root):
    files = sorted((repo_root / "data").glob("workspace_fk_ins94.3_*.npz"))
    if not files:
        return float("nan")
    data = np.load(files[0])
    U = data["U"]
    P = data["P"]
    y = P[:, 1]
    c3 = U[:, 2]
    eps = 1e-3
    signs = np.array([_sign(val, eps) for val in y])
    c3_signs = np.array([_sign(val, eps) for val in c3])
    agree = signs == c3_signs
    return float(np.mean(agree))


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


def _rollout_stats(ext, cfg, li, x_base, u_seq, tip_dyn, use_crm_step):
    residuals = []
    exit_flags = []
    errors = []

    x_t = x_base
    for t in range(u_seq.shape[0]):
        u_t = u_seq[t:t + 1]
        if use_crm_step:
            x_tp1 = crm_step(x_t, u_t, li, cfg)
            x_tp1_ref, eta, exit_flag, residual_norm, u0, tau = ext.crm_step_forward(x_t, u_t, li, cfg)
        else:
            x_tp1, eta, exit_flag, residual_norm, u0, tau = ext.crm_step_forward(x_t, u_t, li, cfg)
        exit_flags.append(int(exit_flag.item()))
        residuals.append(float(residual_norm.item()))

        tip_pred = x_tp1.detach()[0, 24:27].cpu().numpy()
        tip_ref = tip_dyn[t]
        errors.append(float(np.linalg.norm(tip_pred - tip_ref)))

        x_t = x_tp1.detach()

    return {
        "max_tip_error": float(np.max(errors)),
        "mean_tip_error": float(np.mean(errors)),
        "max_residual_norm": float(np.max(residuals)),
        "convergence_rate": float(np.mean(np.array(exit_flags) == 0)),
    }


def main():
    repo_root = Path(__file__).resolve().parents[1]
    out_dir = repo_root / "output_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "stress_report_943.csv"

    ext = load_dyn_ext()
    ext_eq = _load_equilibrium_ext()
    x_base, _ = build_state()

    sign_rate = _workspace_sign_agreement(repo_root)

    rows = []
    for fname in sorted(glob.glob(str(repo_root / "data" / "dyn_fk_*.npz"))):
        data = np.load(fname)
        u_seq = torch.tensor(data["currents"], dtype=torch.float64).view(-1, 1, 3)
        tip_dyn = data["tip_dyn"]

        li = torch.tensor([float(np.array(data["insertion_length"]).item())], dtype=torch.float64)
        dt = float(np.array(data["dt"]).item())
        integration_step = float(np.array(data["integration_step_size"]).item()) if "integration_step_size" in data.files else 0.2
        cfg = _build_cfg(ext, integration_step, dt)

        cfg_eq = _build_equilibrium_cfg(ext_eq, integration_step)
        x_eq = _equilibrium_seed_state(ext_eq, cfg_eq, li, u_seq[0, 0].cpu().numpy())

        modes = {
            "zero": (x_base, False),
            "warm": (x_base, True),
            "workspace": (x_base, False),
            "equilibrium_seed": (x_eq, True),
        }
        for mode, (x0, use_crm_step) in modes.items():
            stats = _rollout_stats(ext, cfg, li, x0, u_seq, tip_dyn, use_crm_step)
            rows.append({
                "dataset_name": Path(fname).name,
                "init_mode": mode,
                "max_tip_error": stats["max_tip_error"],
                "mean_tip_error": stats["mean_tip_error"],
                "max_residual_norm": stats["max_residual_norm"],
                "convergence_rate": stats["convergence_rate"],
                "y_c3_sign_agreement_rate": float("nan"),
            })

    rows.append({
        "dataset_name": "workspace_fk_ins94.3",
        "init_mode": "workspace",
        "max_tip_error": float("nan"),
        "mean_tip_error": float("nan"),
        "max_residual_norm": float("nan"),
        "convergence_rate": float("nan"),
        "y_c3_sign_agreement_rate": sign_rate,
    })

    with open(out_path, "w", encoding="ascii", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "dataset_name",
                "init_mode",
                "max_tip_error",
                "mean_tip_error",
                "max_residual_norm",
                "convergence_rate",
                "y_c3_sign_agreement_rate",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
