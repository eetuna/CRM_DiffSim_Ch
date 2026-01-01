import csv
import json
from datetime import datetime
from pathlib import Path


def main():
    repo_root = Path(__file__).resolve().parents[1]
    csv_path = repo_root / "output_data" / "forward_stability_sweep_captured.csv"
    if not csv_path.exists():
        raise SystemExit(f"missing {csv_path}")

    safe_rows = []
    with csv_path.open("r", encoding="ascii") as f:
        reader = csv.DictReader(f)
        for row in reader:
            nan_detected = row["nan_detected"].lower() == "true"
            unbounded = row["unbounded_detected"].lower() == "true"
            if not nan_detected and not unbounded:
                safe_rows.append(row)

    if not safe_rows:
        raise SystemExit("no SAFE rows found in captured CSV")

    umax_safe = max(float(r["umax"]) for r in safe_rows)
    d_umax_safe = max(float(r["d_umax"]) for r in safe_rows)

    out_dir = repo_root / "docs" / "control"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "safe_bounds.json"

    payload = {
        "umax_safe": umax_safe,
        "d_umax_safe": d_umax_safe,
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "csv_path": str(csv_path),
    }

    with out_path.open("w", encoding="ascii") as f:
        json.dump(payload, f, indent=2)

    print(f"umax_safe={umax_safe} d_umax_safe={d_umax_safe}")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
