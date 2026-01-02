import os
import sys

import numpy as np
import torch

repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if repo_root not in sys.path:
    sys.path.insert(0, repo_root)
python_root = os.path.join(repo_root, "python")
if python_root not in sys.path:
    sys.path.insert(0, python_root)

from crm_diffsims.control.ilqr import tip_position
from crm_diffsims.dynamics.step import build_cfg, build_state, load_dyn_ext
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward


def test_open_loop_replay_matches_dataset():
    torch.manual_seed(0)
    np.random.seed(0)

    data = np.load("data/dyn_fk_ramp_circle1_hold1.npz")
    currents = data["currents"]
    tip_dyn = data["tip_dyn"]
    li_val = float(data.get("insertion_length", 94.3))

    ext = load_dyn_ext()
    cfg_dyn = build_cfg(ext)
    x0, _ = build_state()
    x0 = x0.clone()
    x0[0, 24:27] = torch.tensor(tip_dyn[0], dtype=x0.dtype)
    li = torch.tensor([li_val], dtype=torch.float64)

    n_steps = min(20, currents.shape[0])
    u_seq = torch.tensor(currents[:n_steps], dtype=torch.float64).view(n_steps, 1, 3)
    tip_preds = [tip_position(x0)[0].detach().cpu().numpy()]
    for t in range(n_steps):
        x_tp1, _ = crm_step_v1_3_forward(x0, u_seq[t : t + 1], li, cfg_dyn)
        x0 = x_tp1
        tip_preds.append(tip_position(x_tp1)[0].detach().cpu().numpy())

    tip_preds = np.asarray(tip_preds)
    tip_ref = tip_dyn[: tip_preds.shape[0]]
    errors = np.linalg.norm(tip_preds - tip_ref, axis=1)
    assert float(errors.mean()) < 12.0
    assert float(errors.max()) < 40.0
