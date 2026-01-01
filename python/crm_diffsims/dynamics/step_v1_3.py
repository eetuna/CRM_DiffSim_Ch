"""v1.3 dynamics step scaffolding (implicit backward TODO)."""

import torch

from crm_diffsims.dynamics.step import load_dyn_ext


class CRMStepV13Function(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x_t, u_t, Li, cfg, step_idx=-1):
        ext = load_dyn_ext()
        x_tp1, eta, solver_exit, residual_norm, u0, tau = ext.crm_step_forward(x_t, u_t, Li, cfg)
        ctx.save_for_backward(x_t, u_t, Li, eta, solver_exit, residual_norm, u0, tau)
        ctx.cfg = cfg
        ctx.step_idx = step_idx
        return x_tp1

    @staticmethod
    def backward(ctx, grad_output):
        ext = load_dyn_ext()
        x_t, u_t, Li, eta, solver_exit, residual_norm, u0, tau = ctx.saved_tensors
        res_thresh = ctx.cfg.residual_threshold

        solver_exit_cpu = solver_exit.cpu()
        residual_cpu = residual_norm.cpu()
        for b in range(solver_exit_cpu.numel()):
            exit_flag = int(solver_exit_cpu[b].item())
            res_val = float(residual_cpu[b].item())
            if exit_flag != 0 or res_val > res_thresh:
                raise RuntimeError(
                    "[crm_step_v1_3_backward] non-convergence: batch="
                    f"{b} step={ctx.step_idx} solver_exit_flag={exit_flag} residual_norm={res_val}"
                )

        vjp_x, vjp_u = ext.crm_dyn_vjp_v1_3(grad_output, eta, x_t, u_t, Li, ctx.cfg)
        return vjp_x, vjp_u, None, None, None


def crm_step_v1_3(x_t, u_t, Li, cfg_dyn):
    return CRMStepV13Function.apply(x_t, u_t, Li, cfg_dyn)


def crm_step_v1_3_forward(x_t, u_t, Li, cfg_dyn):
    ext = load_dyn_ext()
    x_tp1, eta, solver_exit, residual_norm, u0, tau = ext.crm_step_forward(x_t, u_t, Li, cfg_dyn)
    cache = {
        "eta": eta.detach(),
        "solver_exit": solver_exit.detach(),
        "residual_norm": residual_norm.detach(),
        "u0": u0.detach(),
        "tau": tau.detach(),
    }
    return x_tp1, cache
