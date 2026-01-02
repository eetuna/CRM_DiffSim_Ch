import os
import sys
import threading
import time

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
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    build_dir = os.path.join(repo_root, "build_torch_ext_dyn")
    os.environ["TORCH_EXTENSIONS_DIR"] = build_dir
    print(f"TORCH_EXTENSIONS_DIR={build_dir}", flush=True)
    lock_path = os.path.join(build_dir, "lock")
    if os.path.exists(lock_path):
        print(f"Lock exists: {lock_path}", flush=True)
        raise SystemExit(2)
    stop_event = threading.Event()

    def _heartbeat():
        start = time.time()
        warned = False
        while not stop_event.is_set():
            elapsed = time.time() - start
            print(f"compiling/ext-loading... elapsed_sec={elapsed:.1f}", flush=True)
            if elapsed > 120.0 and not warned:
                lock_hint = os.path.join(build_dir, "lock")
                print(
                    f"wait >120s; check lock at {lock_hint} or remove ~/.cache/torch_extensions/*",
                    flush=True,
                )
                warned = True
            stop_event.wait(2.0)

    thread = threading.Thread(target=_heartbeat, daemon=True)
    thread.start()
    try:
        ext = load_dyn_ext()
    finally:
        stop_event.set()
        thread.join(timeout=1.0)
    cfg = build_cfg(ext)
    x_t, li = build_state()
    u_t = torch.zeros((1, 1, 3), dtype=torch.float64)
    x_tp1 = crm_step(x_t, u_t, li, cfg)
    if not torch.isfinite(x_tp1).all().item():
        raise RuntimeError("Warmup produced non-finite state")
    print("WARMUP OK")


if __name__ == "__main__":
    main()
