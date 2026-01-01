"""Plot helpers for CRM control demos."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def _ensure_parent(path: str | Path) -> Path:
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    return out_path


def _save_both(fig, out_path: str | Path) -> None:
    out_path = _ensure_parent(out_path)
    fig.savefig(out_path.with_suffix(".png"), dpi=200)
    fig.savefig(out_path.with_suffix(".pdf"))


def save_tip_trajectory_plot(tip_xyz: np.ndarray, target_xyz: np.ndarray, *, title: str, out_path: str | Path) -> None:
    try:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting") from exc

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(target_xyz[:, 0], target_xyz[:, 1], target_xyz[:, 2], label="target")
    ax.plot(tip_xyz[:, 0], tip_xyz[:, 1], tip_xyz[:, 2], label="cem")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend()
    fig.tight_layout()
    _save_both(fig, out_path)


def save_currents_plot(u_seq: np.ndarray, *, title: str, out_path: str | Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting") from exc

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(u_seq[:, 0], label="i1")
    ax.plot(u_seq[:, 1], label="i2")
    ax.plot(u_seq[:, 2], label="i3")
    ax.set_title(title)
    ax.set_xlabel("time")
    ax.set_ylabel("current")
    ax.legend()
    fig.tight_layout()
    _save_both(fig, out_path)


def save_error_plot(errors: np.ndarray, *, title: str, out_path: str | Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting") from exc

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(errors, label="tip error")
    ax.set_title(title)
    ax.set_xlabel("time")
    ax.set_ylabel("error")
    ax.legend()
    fig.tight_layout()
    _save_both(fig, out_path)
