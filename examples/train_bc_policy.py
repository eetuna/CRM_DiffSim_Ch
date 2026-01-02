import argparse
import os
import sys
from datetime import datetime

if __package__ is None:
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    python_root = os.path.join(repo_root, "python")
    if python_root not in sys.path:
        sys.path.insert(0, python_root)

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, Subset

from crm_diffsims.control.run_report import write_run_report
from crm_diffsims.learning.bc_dataset import BCDataset


class MLPPolicy(nn.Module):
    def __init__(self, input_dim, hidden_sizes):
        super().__init__()
        layers = []
        dim = input_dim
        for h in hidden_sizes:
            layers.append(nn.Linear(dim, h))
            layers.append(nn.ReLU())
            dim = h
        layers.append(nn.Linear(dim, 3))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)


def _parse_hidden(text):
    return [int(v) for v in text.split(",") if v.strip()]


def _evaluate(model, loader, device):
    model.eval()
    errs = []
    with torch.no_grad():
        for x_in, u_gt, _ in loader:
            x_in = x_in.to(device)
            u_gt = u_gt.to(device)
            pred = model(x_in)
            err = torch.norm(pred - u_gt, dim=-1)
            errs.append(err.cpu())
    if not errs:
        return 0.0, 0.0
    errs = torch.cat(errs, dim=0)
    return float(errs.mean().item()), float(errs.max().item())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, required=True)
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--hidden", type=str, default="64,64")
    parser.add_argument("--include-tip", action="store_true")
    parser.add_argument("--include-prev-u", action="store_true")
    parser.add_argument("--smooth-weight", type=float, default=0.1)
    parser.add_argument("--val-split", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--checkpoint-every", type=int, default=0)
    parser.add_argument("--output-model", type=str, default=None)
    parser.add_argument("--output-npz", type=str, default=None)
    parser.add_argument("--report-path", type=str, default=None)
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    dataset = BCDataset(
        args.dataset,
        include_tip=args.include_tip,
        include_prev_u=args.include_prev_u,
        dtype=torch.float32,
    )
    print(
        "dataset_counts: "
        f"total={dataset.meta.n_steps_total} kept={dataset.meta.n_steps_kept} "
        f"dropped={dataset.meta.n_steps_dropped} reasons={dataset.meta.drop_reasons_counts}"
    )

    n_total = len(dataset)
    indices = np.arange(n_total)
    rng = np.random.default_rng(args.seed)
    rng.shuffle(indices)
    split = int(n_total * (1.0 - args.val_split))
    train_idx = indices[:split]
    val_idx = indices[split:] if split < n_total else indices[:0]

    train_set = Subset(dataset, train_idx)
    val_set = Subset(dataset, val_idx)

    train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_set, batch_size=args.batch_size, shuffle=False)
    full_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False)

    hidden_sizes = _parse_hidden(args.hidden)
    model = MLPPolicy(dataset.input_dim, hidden_sizes)
    device = torch.device("cpu")
    model.to(device)

    optim = torch.optim.Adam(model.parameters(), lr=args.lr)
    loss_fn = nn.MSELoss()

    history = {"loss": [], "val_loss": []}

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    model_path = args.output_model or os.path.join("output_data", "models", f"bc_policy_{stamp}.pt")
    os.makedirs(os.path.dirname(model_path), exist_ok=True)

    for epoch in range(args.epochs):
        model.train()
        epoch_loss = 0.0
        n_seen = 0
        for x_in, u_gt, prev_u in train_loader:
            x_in = x_in.to(device)
            u_gt = u_gt.to(device)
            prev_u = prev_u.to(device)

            pred = model(x_in)
            mse = loss_fn(pred, u_gt)
            smooth = loss_fn(pred - prev_u, torch.zeros_like(pred))
            loss = mse + args.smooth_weight * smooth

            optim.zero_grad()
            loss.backward()
            optim.step()

            batch = x_in.shape[0]
            epoch_loss += loss.item() * batch
            n_seen += batch

        history["loss"].append(epoch_loss / max(n_seen, 1))

        if len(val_loader) > 0:
            model.eval()
            val_loss = 0.0
            n_val = 0
            with torch.no_grad():
                for x_in, u_gt, prev_u in val_loader:
                    x_in = x_in.to(device)
                    u_gt = u_gt.to(device)
                    prev_u = prev_u.to(device)
                    pred = model(x_in)
                    mse = loss_fn(pred, u_gt)
                    smooth = loss_fn(pred - prev_u, torch.zeros_like(pred))
                    loss = mse + args.smooth_weight * smooth
                    batch = x_in.shape[0]
                    val_loss += loss.item() * batch
                    n_val += batch
            history["val_loss"].append(val_loss / max(n_val, 1))
        else:
            history["val_loss"].append(0.0)

        if args.checkpoint_every > 0 and (epoch + 1) % args.checkpoint_every == 0:
            checkpoint_path = model_path.replace(".pt", f"_epoch{epoch + 1}.pt")
            torch.save(
                {
                    "state_dict": model.state_dict(),
                    "input_dim": dataset.input_dim,
                    "hidden_sizes": hidden_sizes,
                    "include_tip": args.include_tip,
                    "include_prev_u": args.include_prev_u,
                    "umax": dataset.meta.umax,
                    "d_umax": dataset.meta.d_umax,
                },
                checkpoint_path,
            )

    mean_error, max_error = _evaluate(model, full_loader, device)

    torch.save(
        {
            "state_dict": model.state_dict(),
            "input_dim": dataset.input_dim,
            "hidden_sizes": hidden_sizes,
            "include_tip": args.include_tip,
            "include_prev_u": args.include_prev_u,
            "umax": dataset.meta.umax,
            "d_umax": dataset.meta.d_umax,
        },
        model_path,
    )

    out_npz = args.output_npz or os.path.join("output_data", f"bc_train_{stamp}.npz")
    np.savez(
        out_npz,
        loss=np.asarray(history["loss"], dtype=np.float64),
        val_loss=np.asarray(history["val_loss"], dtype=np.float64),
        mean_error=mean_error,
        max_error=max_error,
        dataset=args.dataset,
        model_path=model_path,
    )

    report_path = args.report_path or os.path.join(
        "docs",
        "control",
        "run_reports",
        f"train_bc_policy_{stamp}.md",
    )
    write_run_report(
        report_path=report_path,
        command=" ".join(sys.argv),
        dataset_name=dataset.meta.dataset_name,
        start_idx=dataset.meta.start_idx,
        end_idx=dataset.meta.end_idx,
        n_total=dataset.meta.n_total,
        li_mm=dataset.meta.li,
        dt=dataset.meta.dt,
        umax=dataset.meta.umax,
        d_umax=dataset.meta.d_umax,
        mean_error=mean_error,
        max_error=max_error,
        artifact_npz=out_npz,
        plot_prefix=None,
        init_from_ramp=dataset.meta.ramp_len > 0,
        ramp_len=dataset.meta.ramp_len,
        max_wall_hit=False,
        extra_entries={
            "n_steps_total": dataset.meta.n_steps_total,
            "n_steps_kept": dataset.meta.n_steps_kept,
            "n_steps_dropped": dataset.meta.n_steps_dropped,
            "drop_reasons_counts": dataset.meta.drop_reasons_counts,
        },
    )

    print(f"Saved model to {model_path}")
    print(f"Training metrics: {out_npz}")
    print(f"Run report: {report_path}")


if __name__ == "__main__":
    main()
