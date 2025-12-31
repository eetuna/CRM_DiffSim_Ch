"""Finite-difference autograd bridge for CRM step dynamics."""

from pathlib import Path

import torch
from torch.utils.cpp_extension import load

_EXT = None


def load_dyn_ext():
    global _EXT
    if _EXT is not None:
        return _EXT
    repo_root = Path(__file__).resolve().parents[3]
    build_dir = repo_root / "build_torch_ext_dyn"
    build_dir.mkdir(parents=True, exist_ok=True)
    sources = [
        str(repo_root / "src" / "CRM_TorchDynamics.cpp"),
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
    _EXT = load(
        name="crm_dynamics_ext",
        sources=sources,
        extra_include_paths=include_dirs,
        extra_cflags=["-O3", "-std=c++17"],
        build_directory=str(build_dir),
        verbose=False,
    )
    return _EXT


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


def build_cfg(ext):
    repo_root = Path(__file__).resolve().parents[3]
    cath_params_path = repo_root / "catheterdata" / "CatheterParameterSet_1_dyn.txt"
    cath_config_path = repo_root / "catheterdata" / "CatheterSpatialConfiguration_1.txt"

    params = _parse_param_file(cath_params_path)
    act_inertia = _act_inertia_from_params(params)

    damping = [
        12.1761626666366, 12.1761626666366, 284.429938756989,
        0.0304776127617393, 0.0304776127617393, 0.00502712804532508,
    ]

    cfg = ext.DynamicsConfig(
        str(cath_params_path),
        str(cath_config_path),
        0.2,
        [0.0, 0.0, 0.0],
        act_inertia,
        damping,
        0.05,
        1e-2,
        1e-5,
    )
    return cfg


def build_state():
    xf_pre = [
        -0.458414144062750,
        34.411241976876518,
        70.457561147732264,
        0.999932718178103,
        0.009921777042635,
        -0.006009780134551,
        -0.004651734390922,
        0.817579117734723,
        0.575797488368325,
        0.010626405041486,
        -0.575730791763330,
        0.817570262993625,
        -0.015378744286498,
        0.000001280646594,
        -0.000349413951059,
    ]
    p_l = [-0.248418562587657, 17.707660318406560, 46.752162601547091]
    r_l = [
        0.999919687839427, 0.009924211584043, -0.007882125064742,
        -0.003571217614502, 0.817374079004311, 0.576096225796181,
        0.012159945552960, -0.576021809479719, 0.817343875445250,
    ]
    v_l = [0.0, 0.0, 0.0]
    w_l = [0.0, 0.0, 0.0]
    m_l = [0.0, 0.0, 0.0]
    n_l = [0.0, 0.0, 0.0]

    x_t = torch.tensor(
        v_l + w_l + m_l + n_l + p_l + r_l + xf_pre,
        dtype=torch.float64,
    ).view(1, -1)
    li = torch.tensor([94.3], dtype=torch.float64)
    return x_t, li


def build_stable_state(ext, cfg):
    x_t, li = build_state()
    u_t = torch.tensor([[[0.0, 0.0, 0.0]]], dtype=torch.float64)
    x_tp1, eta, exit_flag, res_norm, u0, tau = ext.crm_step_forward(x_t, u_t, li, cfg)
    if int(exit_flag.item()) != 0:
        raise RuntimeError(f"stable_state solver_exit_flag={int(exit_flag.item())}")
    return x_tp1.detach(), li


class CRMStepFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x_t, u_t, Li, cfg, step_idx=-1):
        ext = load_dyn_ext()
        x_tp1, eta, solver_exit, residual_norm, u0, tau = ext.crm_step_forward(x_t, u_t, Li, cfg)
        ctx.save_for_backward(x_t, u_t, Li, eta, solver_exit, residual_norm)
        ctx.cfg = cfg
        ctx.step_idx = step_idx
        return x_tp1

    @staticmethod
    def backward(ctx, grad_output):
        ext = load_dyn_ext()
        x_t, u_t, Li, eta, solver_exit, residual_norm = ctx.saved_tensors
        res_thresh = ctx.cfg.residual_threshold

        solver_exit_cpu = solver_exit.cpu()
        residual_cpu = residual_norm.cpu()
        for b in range(solver_exit_cpu.numel()):
            exit_flag = int(solver_exit_cpu[b].item())
            res_val = float(residual_cpu[b].item())
            if exit_flag != 0 or res_val > res_thresh:
                raise RuntimeError(
                    "[crm_step_backward] non-convergence: batch="
                    f"{b} step={ctx.step_idx} solver_exit_flag={exit_flag} residual_norm={res_val}"
                )

        grad_x = torch.zeros_like(x_t)
        grad_u = torch.zeros_like(u_t)
        fd_eps = 1e-6
        for b in range(x_t.shape[0]):
            g = grad_output[b]

            for j in range(x_t.shape[1]):
                x_plus = x_t.clone()
                x_minus = x_t.clone()
                x_plus[b, j] += fd_eps
                x_minus[b, j] -= fd_eps
                y_plus = ext.crm_step_forward(x_plus, u_t, Li, ctx.cfg)[0][b]
                y_minus = ext.crm_step_forward(x_minus, u_t, Li, ctx.cfg)[0][b]
                grad_x[b, j] = (y_plus - y_minus).dot(g) / (2.0 * fd_eps)

            for j in range(u_t.shape[2]):
                u_plus = u_t.clone()
                u_minus = u_t.clone()
                u_plus[b, 0, j] += fd_eps
                u_minus[b, 0, j] -= fd_eps
                y_plus = ext.crm_step_forward(x_t, u_plus, Li, ctx.cfg)[0][b]
                y_minus = ext.crm_step_forward(x_t, u_minus, Li, ctx.cfg)[0][b]
                grad_u[b, 0, j] = (y_plus - y_minus).dot(g) / (2.0 * fd_eps)

        return grad_x, grad_u, None, None, None


def crm_step(x_t, u_t, Li, cfg_dyn):
    return CRMStepFunction.apply(x_t, u_t, Li, cfg_dyn)
