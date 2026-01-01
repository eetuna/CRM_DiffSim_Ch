"""Linearization helper using v1.3 implicit VJP/JVP."""

from dataclasses import dataclass
from typing import Callable, Tuple

import torch

from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3


@dataclass
class LinearOperator:
    apply: Callable[[torch.Tensor], torch.Tensor]
    apply_transpose: Callable[[torch.Tensor], torch.Tensor]


def _jvp_step(x_t, u_t, li, cfg_dyn, dx, du):
    def func(x_in, u_in):
        return crm_step_v1_3(x_in, u_in, li, cfg_dyn)

    _, jvp = torch.autograd.functional.jvp(
        func,
        (x_t, u_t),
        (dx, du),
        create_graph=False,
        strict=False,
    )
    return jvp


def linearize_v1_3(
    x_t: torch.Tensor,
    u_t: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    *,
    return_dense: bool = False,
) -> Tuple[LinearOperator, LinearOperator] | Tuple[torch.Tensor, torch.Tensor]:
    """Return linear operators A,B or dense A,B (optional).

    A·dx ≈ f(x+dx,u) - f(x,u)
    B·du ≈ f(x,u+du) - f(x,u)

    Note: valid only for small perturbations (trust region).
    """
    x_t = x_t.detach()
    u_t = u_t.detach()

    def apply_a(dx):
        return _jvp_step(x_t, u_t, li, cfg_dyn, dx, torch.zeros_like(u_t))

    def apply_b(du):
        return _jvp_step(x_t, u_t, li, cfg_dyn, torch.zeros_like(x_t), du)

    def apply_a_t(v):
        vjp_x, _ = torch.autograd.functional.vjp(
            lambda x_in: crm_step_v1_3(x_in, u_t, li, cfg_dyn),
            x_t,
            v,
            create_graph=False,
            strict=True,
        )
        return vjp_x

    def apply_b_t(v):
        vjp_u, = torch.autograd.functional.vjp(
            lambda u_in: crm_step_v1_3(x_t, u_in, li, cfg_dyn),
            u_t,
            v,
            create_graph=False,
            strict=True,
        )
        return vjp_u

    if not return_dense:
        return LinearOperator(apply=apply_a, apply_transpose=apply_a_t), LinearOperator(
            apply=apply_b, apply_transpose=apply_b_t
        )

    n_x = x_t.shape[-1]
    n_u = u_t.shape[-1]
    basis_x = torch.eye(n_x, dtype=x_t.dtype, device=x_t.device)
    basis_u = torch.eye(n_u, dtype=u_t.dtype, device=u_t.device)

    a_cols = []
    for i in range(n_x):
        dx = torch.zeros_like(x_t)
        dx[0, :] = basis_x[i]
        a_cols.append(apply_a(dx)[0])
    a = torch.stack(a_cols, dim=1)

    b_cols = []
    for i in range(n_u):
        du = torch.zeros_like(u_t)
        du[0, 0, :] = basis_u[i]
        b_cols.append(apply_b(du)[0])
    b = torch.stack(b_cols, dim=1)

    return a, b
