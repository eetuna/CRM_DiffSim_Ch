from pathlib import Path

import torch
from torch.utils.cpp_extension import load


def _load_ext():
    repo_root = Path(__file__).resolve().parents[1]
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
    return load(
        name="crm_dynamics_ext",
        sources=sources,
        extra_include_paths=include_dirs,
        extra_cflags=["-O3", "-std=c++17"],
        build_directory=str(build_dir),
        verbose=False,
    )


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


def test_step_forward_deterministic():
    torch.manual_seed(0)
    ext = _load_ext()

    repo_root = Path(__file__).resolve().parents[1]
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
        1e-5,
        1e-5,
    )

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
    u_t = torch.tensor([[[0.0, 0.0, 0.1]]], dtype=torch.float64)
    li = torch.tensor([94.3], dtype=torch.float64)

    out1 = ext.crm_step_forward(x_t, u_t, li, cfg)
    out2 = ext.crm_step_forward(x_t, u_t, li, cfg)

    x_tp1_1, eta_1, exit_1, res_1, u0_1, tau_1 = out1
    x_tp1_2, eta_2, exit_2, res_2, u0_2, tau_2 = out2

    assert x_tp1_1.shape == x_t.shape
    assert eta_1.shape[1] == 6

    assert torch.allclose(x_tp1_1, x_tp1_2, atol=1e-12, rtol=0.0)
    assert torch.allclose(eta_1, eta_2, atol=1e-12, rtol=0.0)
    assert torch.equal(exit_1, exit_2)
    assert torch.allclose(res_1, res_2, atol=1e-12, rtol=0.0)
    assert torch.allclose(u0_1, u0_2, atol=1e-12, rtol=0.0)
    assert torch.allclose(tau_1, tau_2, atol=1e-12, rtol=0.0)
