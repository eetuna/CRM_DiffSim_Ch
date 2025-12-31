"""Step timing profiler and flight recorder for CRM dynamics."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import torch

from crm_diffsims.dynamics.step import load_dyn_ext


@dataclass
class StepProfileResult:
    timings: List[float]
    tips: List[List[float]]
    nan_steps: List[int]
    fail_index: Optional[int]
    total_time: float
    mean_step_time: float
    max_step_time: float
    exit_flags: Optional[List[int]]
    residual_norms: Optional[List[float]]


def _safe_tensor_to_list(x: torch.Tensor) -> List[float]:
    return [float(v) for v in x.detach().cpu().tolist()]


def run_rollout_profile(
    x0: torch.Tensor,
    u_seq: torch.Tensor,
    li: torch.Tensor,
    cfg_dyn,
    *,
    label: str,
    save_dir: str | Path,
) -> StepProfileResult:
    """Run a rollout and record per-step timing/tip/flags.

    Stops early if NaN/Inf appears and writes JSON/NPZ artifacts.
    """
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)

    ext = load_dyn_ext()

    timings: List[float] = []
    tips: List[List[float]] = []
    nan_steps: List[int] = []
    exit_flags: List[int] = []
    residual_norms: List[float] = []
    fail_index: Optional[int] = None

    x_t = x0
    start = time.perf_counter()
    for t in range(u_seq.shape[0]):
        t0 = time.perf_counter()
        x_tp1, _, exit_flag, res_norm, _, _ = ext.crm_step_forward(x_t, u_seq[t : t + 1], li, cfg_dyn)
        dt = time.perf_counter() - t0
        timings.append(dt)

        tip = x_tp1[0, 24:27]
        tips.append(_safe_tensor_to_list(tip))
        exit_flags.append(int(exit_flag.item()))
        residual_norms.append(float(res_norm.item()))

        if not torch.isfinite(x_tp1).all().item():
            nan_steps.append(t)
            fail_index = t
            x_t = x_tp1
            break

        x_t = x_tp1

    total_time = time.perf_counter() - start
    mean_step_time = float(np.mean(timings)) if timings else 0.0
    max_step_time = float(np.max(timings)) if timings else 0.0

    payload = {
        "label": label,
        "timings": timings,
        "tips": tips,
        "nan_steps": nan_steps,
        "fail_index": fail_index,
        "exit_flags": exit_flags,
        "residual_norms": residual_norms,
        "total_time": total_time,
        "mean_step_time": mean_step_time,
        "max_step_time": max_step_time,
    }

    json_path = save_path / f"flight_{label}.json"
    npz_path = save_path / f"flight_{label}.npz"

    with json_path.open("w", encoding="ascii") as f:
        json.dump(payload, f, indent=2)

    np.savez(
        npz_path,
        x0=x0.detach().cpu().numpy(),
        u_seq=u_seq.detach().cpu().numpy(),
        tips=np.asarray(tips, dtype=np.float64),
        timings=np.asarray(timings, dtype=np.float64),
        nan_steps=np.asarray(nan_steps, dtype=np.int32),
        exit_flags=np.asarray(exit_flags, dtype=np.int32),
        residual_norms=np.asarray(residual_norms, dtype=np.float64),
    )

    print(
        "Step profile summary: "
        f"label={label} total_time={total_time:.3f}s mean_step={mean_step_time:.4f}s "
        f"max_step={max_step_time:.4f}s fail_index={fail_index} nan={bool(nan_steps)}"
    )

    return StepProfileResult(
        timings=timings,
        tips=tips,
        nan_steps=nan_steps,
        fail_index=fail_index,
        total_time=total_time,
        mean_step_time=mean_step_time,
        max_step_time=max_step_time,
        exit_flags=exit_flags,
        residual_norms=residual_norms,
    )
