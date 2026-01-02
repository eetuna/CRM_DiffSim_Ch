"""Behavior cloning dataset utilities."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import torch
from torch.utils.data import Dataset


TIP_IDX = slice(24, 27)


@dataclass
class BCDatasetMeta:
    dataset_name: str
    start_idx: int
    end_idx: int
    n_total: int
    n_steps_total: int
    n_steps_kept: int
    n_steps_dropped: int
    drop_reasons_counts: dict
    dt: float
    li: float
    umax: float
    d_umax: float
    ramp_len: int


class BCDataset(Dataset):
    def __init__(
        self,
        npz_path: str,
        *,
        include_tip: bool = True,
        include_prev_u: bool = False,
        dtype: torch.dtype = torch.float32,
    ) -> None:
        payload = np.load(npz_path, allow_pickle=True)
        self.x_t = payload["x_t"].astype(np.float32)
        self.target_xyz = payload["target_xyz"].astype(np.float32)
        self.u_t = payload["u_t"].astype(np.float32)
        self.episode_id = payload["episode_id"].astype(np.int32)
        self.step_id = payload["step_id"].astype(np.int32)

        self.include_tip = include_tip
        self.include_prev_u = include_prev_u
        self.dtype = dtype

        prev_u = np.zeros_like(self.u_t)
        for i in range(1, self.u_t.shape[0]):
            if self.episode_id[i] == self.episode_id[i - 1]:
                prev_u[i] = self.u_t[i - 1]
        self.prev_u = prev_u

        dataset_name = str(payload.get("dataset_name", "unknown"))
        drop_reasons = payload.get("drop_reasons_counts", "{}")
        if isinstance(drop_reasons, np.ndarray):
            drop_reasons = drop_reasons.item()
        if isinstance(drop_reasons, bytes):
            drop_reasons = drop_reasons.decode("ascii")
        try:
            drop_reasons = json.loads(drop_reasons)
        except Exception:
            drop_reasons = {}
        self.meta = BCDatasetMeta(
            dataset_name=dataset_name,
            start_idx=int(payload.get("start_idx", 0)),
            end_idx=int(payload.get("end_idx", self.x_t.shape[0])),
            n_total=int(payload.get("N_total", self.x_t.shape[0])),
            n_steps_total=int(payload.get("n_steps_total", self.x_t.shape[0])),
            n_steps_kept=int(payload.get("n_steps_kept", self.x_t.shape[0])),
            n_steps_dropped=int(payload.get("n_steps_dropped", 0)),
            drop_reasons_counts=drop_reasons,
            dt=float(payload.get("dt", 0.0)),
            li=float(payload.get("Li", 0.0)),
            umax=float(payload.get("umax", 0.1)),
            d_umax=float(payload.get("d_umax", 0.01)),
            ramp_len=int(payload.get("ramp_len", 0)),
        )

    def __len__(self) -> int:
        return self.x_t.shape[0]

    def _build_input(self, idx: int) -> np.ndarray:
        if self.include_tip:
            state = self.x_t[idx, TIP_IDX]
        else:
            state = self.x_t[idx]
        parts = [state, self.target_xyz[idx]]
        if self.include_prev_u:
            parts.append(self.prev_u[idx])
        return np.concatenate(parts, axis=0)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        x_in = torch.tensor(self._build_input(idx), dtype=self.dtype)
        u = torch.tensor(self.u_t[idx], dtype=self.dtype)
        prev_u = torch.tensor(self.prev_u[idx], dtype=self.dtype)
        return x_in, u, prev_u

    @property
    def input_dim(self) -> int:
        base = 3 if self.include_tip else self.x_t.shape[1]
        dim = base + 3
        if self.include_prev_u:
            dim += 3
        return dim
