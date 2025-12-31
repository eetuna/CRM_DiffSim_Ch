import os
from pathlib import Path

import torch
from torch.utils.cpp_extension import load


def _load_ext():
    repo_root = Path(__file__).resolve().parents[1]
    build_dir = repo_root / "build_torch_ext"
    build_dir.mkdir(parents=True, exist_ok=True)
    sources = [
        str(repo_root / "src" / "CRM_TorchEquilibrium.cpp"),
        str(repo_root / "src" / "CRM_IVPSolver.cpp"),
        str(repo_root / "src" / "CRM_BVPSolver.cpp"),
        str(repo_root / "src" / "CRM_IVPJacobian.cpp"),
        str(repo_root / "src" / "CRM_ForwardKinematics.cpp"),
        str(repo_root / "src" / "CRM_CatheterClass.cpp"),
        str(repo_root / "src" / "CRM_SupportFunctions.cpp"),
        str(repo_root / "numerical" / "minpack.cpp"),
        str(repo_root / "numerical" / "minpack_DYN_Defs.cpp"),
        str(repo_root / "src" / "CoilDynamics_Defs.cpp"),
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


class CRMEquilibriumFunction(torch.autograd.Function):
    @staticmethod
    def forward(ctx, u, Li, cfg):
        p_tip, deltau0, ftip, localmin, residual_norm = ext.crm_equilibrium_forward(u, Li, cfg)
        ctx.save_for_backward(u, Li, deltau0, ftip, localmin, residual_norm)
        ctx.cfg = cfg
        return p_tip

    @staticmethod
    def backward(ctx, grad_output):
        u, Li, deltau0, ftip, localmin, residual_norm = ctx.saved_tensors
        grad_u = ext.crm_equilibrium_backward(
            grad_output, u, Li, deltau0, ftip, localmin, residual_norm, ctx.cfg
        )
        return grad_u, None, None


def test_tip_position_gradcheck():
    torch.manual_seed(0)
    global ext
    ext = _load_ext()

    repo_root = Path(__file__).resolve().parents[1]
    cath_params_path = str(repo_root / "catheterdata" / "CatheterParameterSet_1_new.txt")
    cath_config_path = str(repo_root / "catheterdata" / "CatheterSpatialConfiguration_1.txt")

    cfg = ext.EquilibriumConfig(
        cath_params_path,
        cath_config_path,
        0.2,
        [0.0, 0.0, 0.0],
        1e-5,
    )

    Li = torch.tensor([94.0], dtype=torch.float64, requires_grad=False)
    points = [
        [0.0, 0.0, 0.0],
        [0.02, -0.01, 0.03],
        [-0.015, 0.01, 0.02],
    ]

    for p in points:
        u = torch.tensor(p, dtype=torch.float64, requires_grad=True).view(1, 1, 3)

        def func(u_in):
            return CRMEquilibriumFunction.apply(u_in, Li, cfg)

        assert torch.autograd.gradcheck(func, (u,), eps=1e-6, atol=1e-4, rtol=1e-3)
