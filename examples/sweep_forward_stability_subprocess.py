import argparse
import csv
import os
import subprocess
import sys


def _parse_list(text):
    return [float(item) for item in text.split(",") if item.strip()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--umax_list", type=str, required=True)
    parser.add_argument("--d_umax_list", type=str, required=True)
    parser.add_argument("--T", type=int, default=50)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    umax_list = _parse_list(args.umax_list)
    d_umax_list = _parse_list(args.d_umax_list)

    out_dir = "output_data"
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, "forward_stability_sweep_captured.csv")

    fieldnames = [
        "umax",
        "d_umax",
        "T",
        "nan_detected",
        "unbounded_detected",
        "return_code",
    ]

    rows = []
    for umax in umax_list:
        for d_umax in d_umax_list:
            cmd = [
                sys.executable,
                "examples/sweep_forward_stability.py",
                "--single",
                "--umax",
                str(umax),
                "--d_umax",
                str(d_umax),
                "--T",
                str(args.T),
                "--seed",
                str(args.seed),
                "--umax_list",
                str(umax),
                "--d_umax_list",
                str(d_umax),
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            combined = (proc.stdout or "") + (proc.stderr or "")
            unbounded = "Coil integration Unbounded!!" in combined
            nan_detected = proc.returncode != 0

            rows.append(
                {
                    "umax": umax,
                    "d_umax": d_umax,
                    "T": args.T,
                    "nan_detected": nan_detected,
                    "unbounded_detected": unbounded,
                    "return_code": proc.returncode,
                }
            )
            print(
                f"captured umax={umax:.3f} d_umax={d_umax:.3f} "
                f"nan={nan_detected} unbounded={unbounded}"
            )

    with open(csv_path, "w", newline="", encoding="ascii") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved captured sweep results to {csv_path}")


if __name__ == "__main__":
    main()
