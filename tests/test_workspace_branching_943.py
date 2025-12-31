from pathlib import Path

import numpy as np


def _sign(x, eps):
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def test_workspace_branching_943():
    repo_root = Path(__file__).resolve().parents[1]
    files = sorted((repo_root / "data").glob("workspace_fk_ins94.3_*.npz"))
    assert files, "No workspace_fk_ins94.3_*.npz files found in ./data"

    data = np.load(files[0])
    U = data["U"]
    P = data["P"]
    conv = data["conv"]

    assert U.shape[0] == P.shape[0] == conv.shape[0]
    assert np.isfinite(U).all()
    assert np.isfinite(conv).all()

    mask = conv.astype(bool)
    assert mask.any(), "No converged workspace points found."
    assert np.isfinite(P[mask]).all()

    y = P[mask, 1]
    c3 = U[mask, 2]

    eps = 1e-3
    signs = np.array([_sign(val, eps) for val in y])
    c3_signs = np.array([_sign(val, eps) for val in c3])

    mask_pos = signs > 0
    mask_neg = signs < 0
    mask_zero = signs == 0

    agree = signs == c3_signs
    agree_rate = float(np.mean(agree))

    pos_rate = float(np.mean(agree[mask_pos])) if np.any(mask_pos) else float("nan")
    neg_rate = float(np.mean(agree[mask_neg])) if np.any(mask_neg) else float("nan")
    zero_rate = float(np.mean(agree[mask_zero])) if np.any(mask_zero) else float("nan")

    print(f"dataset={files[0].name}")
    print(f"sign_agreement_rate={agree_rate:.4f}")
    print(f"positive_y_rate={pos_rate:.4f}")
    print(f"negative_y_rate={neg_rate:.4f}")
    print(f"near_zero_y_rate={zero_rate:.4f}")

    ambiguous = float(np.mean(mask_zero))
    print(f"ambiguous_fraction={ambiguous:.4f}")

    assert np.all(np.isfinite([agree_rate, pos_rate, neg_rate, ambiguous]))
