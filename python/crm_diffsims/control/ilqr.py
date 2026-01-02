"""Minimal iLQR/MPC utilities for CRM dynamics."""

from dataclasses import dataclass
import json
import os
import time
from typing import Dict, Iterable, List, Tuple

import torch

from crm_diffsims.dynamics.step import crm_step
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


TIP_POS_IDX = torch.tensor([24, 25, 26])


@dataclass
class ILQRConfig:
    horizon: int = 30
    max_iter: int = 10
    tol_cost: float = 1e-5
    reg: float = 1e-4
    reg_min: float = 1e-6
    reg_max: float = 1e2
    w_tip: float = 1.0
    w_tip_terminal: float = 10.0
    w_u: float = 1e-2
    w_du: float = 1e-1
    max_du: float = 0.01
    max_sign_flip: float = 0.005
    max_u: float = 0.05
    cost_increase_ratio: float = 1.2
    max_line_search_tries: int = 3
    line_search_alphas: Tuple[float, ...] = (1.0, 0.5, 0.25, 0.1)
    linearization_backend: str = "autograd"
    linearization_dense: bool = True
    use_safe_bounds: bool = True


@dataclass
class ILQRResult:
    x_seq: torch.Tensor
    u_seq: torch.Tensor
    cost: float
    converged: bool
    iterations: int
    clamp_hits: int
    sign_flip_hits: int
    clamp_rate: float
    sign_flip_rate: float
    nan_detected: bool
    first_iter_u_seq: torch.Tensor | None
    timed_out: bool


def tip_position(x_t: torch.Tensor) -> torch.Tensor:
    return x_t[:, TIP_POS_IDX]


def _maybe_apply_safe_bounds(cfg: ILQRConfig) -> None:
    if not cfg.use_safe_bounds:
        return
    path = "docs/control/safe_bounds.json"
    if not os.path.exists(path):
        return
    try:
        with open(path, "r", encoding="ascii") as f:
            bounds = json.load(f)
    except Exception:
        return
    umax = float(bounds.get("umax_safe", cfg.max_u))
    d_umax = float(bounds.get("d_umax_safe", cfg.max_du))
    cfg.max_u = min(cfg.max_u, umax)
    cfg.max_du = min(cfg.max_du, d_umax)


def _step_forward(
    x_t: torch.Tensor, u_t: torch.Tensor, li: torch.Tensor, cfg_dyn, backend: str
) -> Tuple[torch.Tensor, bool]:
    if backend == "v1_3":
        x_tp1, cache = crm_step_v1_3_forward(x_t, u_t, li, cfg_dyn)
        solver_exit = int(cache["solver_exit"].item())
        residual = float(cache["residual_norm"].item())
        ok = solver_exit == 0 and residual <= cfg_dyn.residual_threshold
        return x_tp1, ok
    return crm_step(x_t, u_t, li, cfg_dyn), True


def rollout(
    x0: torch.Tensor, u_seq: torch.Tensor, li: torch.Tensor, cfg_dyn, backend: str
) -> Tuple[torch.Tensor, torch.Tensor, bool]:
    x_seq = [x0]
    tip_seq = [tip_position(x0)]
    x_t = x0
    ok = True
    for t in range(u_seq.shape[0]):
        x_t, step_ok = _step_forward(x_t, u_seq[t : t + 1], li, cfg_dyn, backend)
        if not step_ok or not torch.isfinite(x_t).all().item():
            ok = False
            break
        x_seq.append(x_t)
        tip_seq.append(tip_position(x_t))
    return torch.cat(x_seq, dim=0), torch.cat(tip_seq, dim=0), ok


def _running_cost(tip: torch.Tensor, target: torch.Tensor, u_t: torch.Tensor) -> torch.Tensor:
    return (tip - target).pow(2).sum(dim=-1) + u_t.pow(2).sum(dim=-1)


def _terminal_cost(tip: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return (tip - target).pow(2).sum(dim=-1)


def total_cost(
    tip_seq: torch.Tensor,
    u_seq: torch.Tensor,
    targets: torch.Tensor,
    cfg: ILQRConfig,
) -> torch.Tensor:
    cost = torch.zeros((), dtype=tip_seq.dtype, device=tip_seq.device)
    for t in range(u_seq.shape[0]):
        tip = tip_seq[t + 1]
        target = targets[t]
        u_t = u_seq[t, 0]
        running = cfg.w_tip * _running_cost(tip, target, u_t)
        cost = cost + running
        if t > 0:
            du = u_seq[t, 0] - u_seq[t - 1, 0]
            cost = cost + cfg.w_du * du.pow(2).sum(dim=-1)
    terminal = cfg.w_tip_terminal * _terminal_cost(tip_seq[-1], targets[-1])
    return cost + terminal


def _linearize(
    x_t: torch.Tensor,
    u_t: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    backend: str,
    dense: bool,
) -> Tuple[torch.Tensor, torch.Tensor]:
    if backend == "v1_3":
        from crm_diffsims.dynamics.linearize_v1_3 import linearize_v1_3

        if not dense:
            raise RuntimeError("iLQR requires dense A,B; set linearization_dense=True.")
        return linearize_v1_3(x_t, u_t, li, cfg_dyn, return_dense=True)

    x_t = x_t.detach().requires_grad_(True)
    u_t = u_t.detach().requires_grad_(True)

    def _dyn_x(x_in):
        return crm_step(x_in, u_t, li, cfg_dyn)

    def _dyn_u(u_in):
        return crm_step(x_t, u_in, li, cfg_dyn)

    a = torch.autograd.functional.jacobian(_dyn_x, x_t, create_graph=False)
    b = torch.autograd.functional.jacobian(_dyn_u, u_t, create_graph=False)

    a = a[0, :, 0, :]
    b = b[0, :, 0, 0, :]
    return a, b


def _cost_derivatives(
    x_t: torch.Tensor,
    u_t: torch.Tensor,
    target: torch.Tensor,
    u_prev: torch.Tensor,
    u_next: torch.Tensor,
    cfg: ILQRConfig,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    n_x = x_t.shape[-1]
    n_u = u_t.shape[-1]
    l_x = torch.zeros(n_x, dtype=x_t.dtype, device=x_t.device)
    l_xx = torch.zeros((n_x, n_x), dtype=x_t.dtype, device=x_t.device)
    l_u = 2.0 * cfg.w_u * u_t
    l_uu = 2.0 * cfg.w_u * torch.eye(n_u, dtype=x_t.dtype, device=x_t.device)

    tip = tip_position(x_t.unsqueeze(0))[0]
    tip_err = tip - target
    l_x[TIP_POS_IDX] = 2.0 * cfg.w_tip * tip_err
    for idx in TIP_POS_IDX.tolist():
        l_xx[idx, idx] = 2.0 * cfg.w_tip

    if u_prev is not None:
        du = u_t - u_prev
        l_u = l_u + 2.0 * cfg.w_du * du
        l_uu = l_uu + 2.0 * cfg.w_du * torch.eye(n_u, dtype=x_t.dtype, device=x_t.device)

    if u_next is not None:
        du_next = u_t - u_next
        l_u = l_u + 2.0 * cfg.w_du * du_next
        l_uu = l_uu + 2.0 * cfg.w_du * torch.eye(n_u, dtype=x_t.dtype, device=x_t.device)

    return l_x, l_u, l_xx, l_uu


def _terminal_derivatives(x_t: torch.Tensor, target: torch.Tensor, cfg: ILQRConfig) -> Tuple[torch.Tensor, torch.Tensor]:
    n_x = x_t.shape[-1]
    v_x = torch.zeros(n_x, dtype=x_t.dtype, device=x_t.device)
    v_xx = torch.zeros((n_x, n_x), dtype=x_t.dtype, device=x_t.device)

    tip = tip_position(x_t.unsqueeze(0))[0]
    tip_err = tip - target
    v_x[TIP_POS_IDX] = 2.0 * cfg.w_tip_terminal * tip_err
    for idx in TIP_POS_IDX.tolist():
        v_xx[idx, idx] = 2.0 * cfg.w_tip_terminal
    return v_x, v_xx


def _apply_constraints(u_prev: torch.Tensor, u_t: torch.Tensor, cfg: ILQRConfig) -> torch.Tensor:
    delta = u_t - u_prev
    delta_clamped = torch.clamp(delta, -cfg.max_du, cfg.max_du)
    clamp_hit = torch.any(delta_clamped != delta).item()
    u_next = u_prev + delta_clamped
    sign_flip = (u_prev * u_next) < 0.0
    sign_flip_hit = sign_flip.any().item()
    if sign_flip_hit:
        delta_clamped = torch.where(
            sign_flip,
            torch.clamp(delta_clamped, -cfg.max_sign_flip, cfg.max_sign_flip),
            delta_clamped,
        )
        u_next = u_prev + delta_clamped
    u_next = torch.clamp(u_next, -cfg.max_u, cfg.max_u)
    return u_next, int(clamp_hit), int(sign_flip_hit)


def tip_index_sanity_check(
    x0: torch.Tensor,
    u0: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    epsilon: float = 1e-3,
    backend: str = "autograd",
) -> Dict[str, float]:
    x_next, _ = _step_forward(x0, u0, li, cfg_dyn, backend)
    u_pert = u0.clone()
    u_pert[0, 0, 0] += epsilon
    x_next_pert, _ = _step_forward(x0, u_pert, li, cfg_dyn, backend)

    tip = tip_position(x_next)[0]
    tip_pert = tip_position(x_next_pert)[0]
    delta = tip_pert - tip

    stats = {
        "tip_finite": float(torch.isfinite(tip).all().item()),
        "tip_delta_norm": float(delta.norm().item()),
        "tip_norm": float(tip.norm().item()),
    }

    print(
        "tip_index_sanity_check: "
        f"tip={tip.tolist()} tip_delta_norm={stats['tip_delta_norm']:.3e} "
        f"finite={bool(stats['tip_finite'])}"
    )
    return stats


def ilqr_solve(
    x0: torch.Tensor,
    u_init: torch.Tensor,
    targets: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    cfg: ILQRConfig,
    max_wall_sec: float | None = None,
    start_time: float | None = None,
) -> ILQRResult:
    _maybe_apply_safe_bounds(cfg)
    u_seq = u_init.clone()
    reg = cfg.reg

    x_seq, tip_seq, rollout_ok = rollout(x0, u_seq, li, cfg_dyn, cfg.linearization_backend)
    cost = total_cost(tip_seq, u_seq, targets, cfg)
    clamp_hits = 0
    sign_flip_hits = 0
    nan_detected = not torch.isfinite(cost).item() or not rollout_ok
    first_iter_u_seq = None
    timed_out = False

    def _check_timeout() -> bool:
        nonlocal timed_out
        if max_wall_sec is not None and start_time is not None:
            if time.time() - start_time > max_wall_sec:
                timed_out = True
                return True
        return False

    converged = False
    iter_clamp_hits = 0
    iter_sign_hits = 0
    for it in range(cfg.max_iter):
        if _check_timeout():
            break
        n_x = x_seq.shape[-1]
        n_u = u_seq.shape[-1]
        k_seq: List[torch.Tensor] = []
        kff_seq: List[torch.Tensor] = []

        v_x, v_xx = _terminal_derivatives(x_seq[-1], targets[-1], cfg)

        for t in reversed(range(u_seq.shape[0])):
            if _check_timeout():
                break
            a_t, b_t = _linearize(
                x_seq[t : t + 1],
                u_seq[t : t + 1],
                li,
                cfg_dyn,
                cfg.linearization_backend,
                cfg.linearization_dense,
            )
            u_prev = u_seq[t - 1, 0] if t > 0 else None
            u_next = u_seq[t + 1, 0] if t < u_seq.shape[0] - 1 else None
            l_x, l_u, l_xx, l_uu = _cost_derivatives(
                x_seq[t],
                u_seq[t, 0],
                targets[t],
                u_prev,
                u_next,
                cfg,
            )

            q_x = l_x + a_t.T @ v_x
            q_u = l_u + b_t.T @ v_x
            q_xx = l_xx + a_t.T @ v_xx @ a_t
            q_ux = b_t.T @ v_xx @ a_t
            q_uu = l_uu + b_t.T @ v_xx @ b_t
            q_uu = q_uu + reg * torch.eye(n_u, dtype=q_uu.dtype, device=q_uu.device)

            q_uu_inv = torch.linalg.inv(q_uu)
            k = -q_uu_inv @ q_ux
            kff = -q_uu_inv @ q_u

            v_x = q_x + k.T @ q_uu @ kff + k.T @ q_u + q_ux.T @ kff
            v_xx = q_xx + k.T @ q_uu @ k + k.T @ q_ux + q_ux.T @ k
            v_xx = 0.5 * (v_xx + v_xx.T)

            k_seq.insert(0, k)
            kff_seq.insert(0, kff)

        if timed_out:
            break

        improved = False
        prev_cost = float(cost)
        iter_clamp_hits = 0
        iter_sign_hits = 0
        for _ in range(cfg.max_line_search_tries):
            for alpha in cfg.line_search_alphas:
                if _check_timeout():
                    break
                x_new = x0
                u_new = torch.zeros_like(u_seq)
                clamp_hits = 0
                sign_flip_hits = 0
                nan_rollout = False
                for t in range(u_seq.shape[0]):
                    if _check_timeout():
                        break
                    dx = (x_new - x_seq[t]).squeeze(0)
                    du = kff_seq[t] * alpha + k_seq[t] @ dx
                    u_prop = (u_seq[t, 0] + du).unsqueeze(0)
                    if t > 0:
                        u_prop, clamp_hit, sign_hit = _apply_constraints(u_new[t - 1, 0], u_prop[0], cfg)
                        clamp_hits += clamp_hit
                        sign_flip_hits += sign_hit
                        u_prop = u_prop.unsqueeze(0)
                    u_new[t, 0] = u_prop[0]
                    x_new, step_ok = _step_forward(x_new, u_new[t : t + 1], li, cfg_dyn, cfg.linearization_backend)
                    if not step_ok or not torch.isfinite(x_new).all().item():
                        nan_rollout = True
                        break

                if nan_rollout or timed_out:
                    nan_detected = True
                    continue

                x_cand, tip_cand, cand_ok = rollout(x0, u_new, li, cfg_dyn, cfg.linearization_backend)
                cost_cand = total_cost(tip_cand, u_new, targets, cfg)
                if not cand_ok or not torch.isfinite(cost_cand).item():
                    nan_detected = True
                    continue

                if cost_cand < cost and cost_cand <= cost * cfg.cost_increase_ratio:
                    if it == 0 and first_iter_u_seq is None:
                        first_iter_u_seq = u_new.detach().clone()
                    u_seq = u_new
                    x_seq = x_cand
                    tip_seq = tip_cand
                    cost = cost_cand
                    iter_clamp_hits = clamp_hits
                    iter_sign_hits = sign_flip_hits
                    improved = True
                    break
                if cost_cand > cost * cfg.cost_increase_ratio:
                    break
            if improved or timed_out:
                break

        if timed_out:
            break

        if not improved:
            reg = min(reg * 10.0, cfg.reg_max)
        else:
            reg = max(reg * 0.5, cfg.reg_min)

        clamp_rate = iter_clamp_hits / max(1, u_seq.shape[0] - 1)
        sign_rate = iter_sign_hits / max(1, u_seq.shape[0] - 1)
        print(
            f"iLQR iter {it + 1}: cost={float(cost):.6e} reg={reg:.3e} "
            f"clamp_rate={clamp_rate:.2%} sign_flip_rate={sign_rate:.2%}"
        )
        if improved and abs(prev_cost - float(cost)) <= cfg.tol_cost:
            converged = True
            break

    clamp_rate = iter_clamp_hits / max(1, u_seq.shape[0] - 1)
    sign_rate = iter_sign_hits / max(1, u_seq.shape[0] - 1)
    return ILQRResult(
        x_seq=x_seq,
        u_seq=u_seq,
        cost=float(cost),
        converged=converged,
        iterations=it + 1,
        clamp_hits=iter_clamp_hits,
        sign_flip_hits=iter_sign_hits,
        clamp_rate=clamp_rate,
        sign_flip_rate=sign_rate,
        nan_detected=nan_detected,
        first_iter_u_seq=first_iter_u_seq,
        timed_out=timed_out,
    )


def run_mpc(
    x0: torch.Tensor,
    targets: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    cfg: ILQRConfig,
    u_warm: torch.Tensor,
    max_wall_sec: float | None = None,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, List[Dict[str, float]]]:
    _maybe_apply_safe_bounds(cfg)
    n_steps = targets.shape[0] - 1
    u_seq = u_warm.clone()
    x_roll = [x0]
    u_roll = []
    stats = []
    start_time = time.time()
    for t in range(n_steps):
        horizon_targets = targets[t : t + cfg.horizon + 1]
        if horizon_targets.shape[0] < cfg.horizon + 1:
            pad = horizon_targets[-1:].repeat(cfg.horizon + 1 - horizon_targets.shape[0], 1)
            horizon_targets = torch.cat([horizon_targets, pad], dim=0)
        result = ilqr_solve(
            x_roll[-1],
            u_seq,
            horizon_targets,
            li,
            cfg_dyn,
            cfg,
            max_wall_sec=max_wall_sec,
            start_time=start_time,
        )
        stats.append(
            {
                "clamp_hits": result.clamp_hits,
                "sign_flip_hits": result.sign_flip_hits,
                "clamp_rate": result.clamp_rate,
                "sign_flip_rate": result.sign_flip_rate,
                "nan_detected": result.nan_detected,
                "first_iter_u_seq": result.first_iter_u_seq,
                "timed_out": result.timed_out,
                "cost": result.cost,
                "converged": result.converged,
            }
        )
        print(
            f"MPC step {t + 1}/{n_steps}: "
            f"clamp_rate={result.clamp_rate:.2%} sign_flip_rate={result.sign_flip_rate:.2%}"
        )
        u_apply = result.u_seq[0:1]
        u_roll.append(u_apply)
        x_next, step_ok = _step_forward(x_roll[-1], u_apply, li, cfg_dyn, cfg.linearization_backend)
        if not step_ok or not torch.isfinite(x_next).all().item():
            break
        x_roll.append(x_next)
        u_seq = torch.cat([result.u_seq[1:], result.u_seq[-1:]], dim=0)
        if result.timed_out:
            break
        if result.nan_detected:
            break
        if max_wall_sec is not None and (time.time() - start_time) > max_wall_sec:
            break
    return torch.cat(x_roll, dim=0), torch.cat(u_roll, dim=0), u_seq, stats
