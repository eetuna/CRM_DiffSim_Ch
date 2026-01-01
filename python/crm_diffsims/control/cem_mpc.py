"""Jacobian-free CEM MPC controller using CRM forward rollouts only."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch

from crm_diffsims.dynamics.step import crm_step


TIP_POS_IDX = torch.tensor([24, 25, 26])


@dataclass
class CEMConfig:
    horizon: int = 20
    mpc_steps: int = 50
    num_samples: int = 64
    num_elites: int = 8
    cem_iters: int = 3
    u_std_init: float = 0.02
    w_tip: float = 1.0
    w_tip_terminal: float = 10.0
    w_u: float = 1e-2
    w_du: float = 1e-1
    umax_safe: float = 0.1
    d_umax_safe: float = 0.01
    max_sign_flip: float = 0.005
    seed: int = 0


@dataclass
class CEMResult:
    x_roll: torch.Tensor
    u_roll: torch.Tensor
    stats: List[Dict[str, float]]


def _load_safe_bounds() -> Tuple[float, float]:
    repo_root = Path(__file__).resolve().parents[3]
    safe_path = repo_root / "docs" / "control" / "safe_bounds.json"
    if not safe_path.exists():
        return 0.1, 0.01
    with safe_path.open("r", encoding="ascii") as f:
        payload = json.load(f)
    return float(payload.get("umax_safe", 0.1)), float(payload.get("d_umax_safe", 0.01))


def tip_position(x_t: torch.Tensor) -> torch.Tensor:
    return x_t[:, TIP_POS_IDX]


def _clamp_sequence(u_seq: torch.Tensor, umax: float, d_umax: float, max_sign_flip: float) -> Tuple[torch.Tensor, int, int]:
    u_clamped = u_seq.clone()
    u_hits = 0
    du_hits = 0
    for t in range(u_clamped.shape[0]):
        before = u_clamped[t].clone()
        u_clamped[t] = torch.clamp(u_clamped[t], -umax, umax)
        if not torch.allclose(before, u_clamped[t]):
            u_hits += 1
        if t > 0:
            delta = u_clamped[t] - u_clamped[t - 1]
            delta_clamped = torch.clamp(delta, -d_umax, d_umax)
            if not torch.allclose(delta, delta_clamped):
                du_hits += 1
            u_next = u_clamped[t - 1] + delta_clamped
            sign_flip = (u_clamped[t - 1] * u_next) < 0.0
            if sign_flip.any():
                delta_clamped = torch.where(
                    sign_flip,
                    torch.clamp(delta_clamped, -max_sign_flip, max_sign_flip),
                    delta_clamped,
                )
                u_next = u_clamped[t - 1] + delta_clamped
            u_clamped[t] = u_next
    return u_clamped, u_hits, du_hits


def _rollout_cost(
    x0: torch.Tensor,
    u_seq: torch.Tensor,
    targets: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    cfg: CEMConfig,
) -> Tuple[torch.Tensor, torch.Tensor, bool]:
    x_t = x0
    tip_seq = []
    cost = torch.zeros((), dtype=x0.dtype)
    nan_detected = False
    for t in range(u_seq.shape[0]):
        x_t = crm_step(x_t, u_seq[t : t + 1], li, cfg_dyn)
        if not torch.isfinite(x_t).all().item():
            nan_detected = True
            cost = cost + 1e6
            break
        tip = tip_position(x_t)[0]
        tip_seq.append(tip)
        target = targets[t]
        cost = cost + cfg.w_tip * (tip - target).pow(2).sum()
        cost = cost + cfg.w_u * u_seq[t, 0].pow(2).sum()
        if t > 0:
            du = u_seq[t, 0] - u_seq[t - 1, 0]
            cost = cost + cfg.w_du * du.pow(2).sum()
    if tip_seq:
        terminal_target = targets[min(len(tip_seq), targets.shape[0] - 1)]
        cost = cost + cfg.w_tip_terminal * (tip_seq[-1] - terminal_target).pow(2).sum()
        tip_seq = torch.stack(tip_seq, dim=0)
    else:
        tip_seq = torch.zeros((0, 3), dtype=x0.dtype)
    return cost, tip_seq, nan_detected


def cem_mpc(
    x0: torch.Tensor,
    targets: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    cfg: CEMConfig,
) -> CEMResult:
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    umax_safe, d_umax_safe = _load_safe_bounds()
    cfg.umax_safe = umax_safe
    cfg.d_umax_safe = d_umax_safe

    mean = torch.zeros((cfg.horizon, 1, 3), dtype=x0.dtype)
    std = torch.full_like(mean, cfg.u_std_init)

    x_roll = [x0]
    u_roll = []
    stats = []

    for step in range(cfg.mpc_steps):
        horizon_targets = targets[step : step + cfg.horizon]
        if horizon_targets.shape[0] < cfg.horizon:
            pad = horizon_targets[-1:].repeat(cfg.horizon - horizon_targets.shape[0], 1)
            horizon_targets = torch.cat([horizon_targets, pad], dim=0)

        all_costs = []
        all_sequences = []
        total_u_hits = 0
        total_du_hits = 0
        for _ in range(cfg.cem_iters):
            samples = mean + std * torch.randn((cfg.num_samples, cfg.horizon, 1, 3), dtype=x0.dtype)
            clamped = []
            u_hits = 0
            du_hits = 0
            for i in range(cfg.num_samples):
                seq, hits_u, hits_du = _clamp_sequence(
                    samples[i],
                    cfg.umax_safe,
                    cfg.d_umax_safe,
                    cfg.max_sign_flip,
                )
                clamped.append(seq)
                u_hits += hits_u
                du_hits += hits_du
            total_u_hits += u_hits
            total_du_hits += du_hits
            clamped = torch.stack(clamped, dim=0)

            costs = []
            for i in range(cfg.num_samples):
                cost, _, nan_flag = _rollout_cost(x_roll[-1], clamped[i], horizon_targets, li, cfg_dyn, cfg)
                if nan_flag:
                    cost = cost + 1e6
                costs.append(cost)
            costs = torch.stack(costs)

            elite_idx = torch.topk(costs, k=cfg.num_elites, largest=False).indices
            elites = clamped[elite_idx]
            mean = elites.mean(dim=0)
            std = elites.std(dim=0).clamp(min=1e-3)

            all_costs = costs
            all_sequences = clamped

        chosen = mean.clone()
        chosen, _, _ = _clamp_sequence(chosen, cfg.umax_safe, cfg.d_umax_safe, cfg.max_sign_flip)
        u0 = chosen[0:1]
        u_roll.append(u0)
        x_next = crm_step(x_roll[-1], u0, li, cfg_dyn)
        x_roll.append(x_next)

        tip_now = tip_position(x_next)[0]
        target_now = horizon_targets[0]
        nan_in_rollout = not torch.isfinite(tip_now).all().item()
        if nan_in_rollout:
            final_err = float("nan")
            mean_err = float("nan")
        else:
            final_err = float(torch.norm(tip_now - target_now).item())
            mean_err = float(torch.norm(tip_position(x_roll[-1])[0] - target_now).item())

        total_candidates = float(cfg.num_samples * cfg.cem_iters * cfg.horizon)
        clamp_rate_u = total_u_hits / max(total_candidates, 1.0)
        clamp_rate_du = total_du_hits / max(total_candidates, 1.0)
        stats.append(
            {
                "mean_tip_error": mean_err,
                "final_tip_error": final_err,
                "clamp_rate_u": clamp_rate_u,
                "clamp_rate_du": clamp_rate_du,
                "u0": float(u0[0, 0, 0].item()),
                "nan_detected": nan_in_rollout,
            }
        )

        mean = torch.cat([mean[1:], mean[-1:]], dim=0)
        if nan_in_rollout:
            break

    return CEMResult(x_roll=torch.cat(x_roll, dim=0), u_roll=torch.cat(u_roll, dim=0), stats=stats)
