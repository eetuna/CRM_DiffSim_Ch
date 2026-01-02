import os
import sys

import torch

if __package__ is None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    python_root = os.path.join(repo_root, "python")
    if python_root not in sys.path:
        sys.path.insert(0, python_root)

from crm_diffsims.dynamics.step import build_cfg, build_state, crm_step, load_dyn_ext


def main():
    torch.manual_seed(0)
    ext = load_dyn_ext()
    cfg = build_cfg(ext)
    x_t, li = build_state()
    u_t = torch.zeros((1, 1, 3), dtype=torch.float64)
    x_tp1 = crm_step(x_t, u_t, li, cfg)
    if not torch.isfinite(x_tp1).all().item():
        raise RuntimeError("Warmup produced non-finite state")
    print("WARMUP OK")


if __name__ == "__main__":
    main()
