# Phase B Closeout

## What Phase B delivered
- Dataset ramp initialization (`--init-from-ramp`) to align the initial state with the dataset regime.
- Planned vs executed separation in artifacts (`executed_*` plus planned tips per MPC step).
- Plot throttling (`--plot-every`) to reduce overhead while still producing PNG+PDF outputs.
- Full-dataset acceptance coverage for circle and lemniscate (baseline replay-only vs controller with ramp init).

## How to reproduce
Circle (baseline replay-only subset):
```
python3 examples/run_ilqr_circle.py --backend v1_3 --full-dataset --mode replay_only_full --start-idx 0 --end-idx 150 --max-wall-sec 300
```

Circle (controller with ramp init):
```
python3 examples/run_ilqr_circle.py --backend v1_3 --full-dataset --mode replay_dataset_mode_full --start-idx 0 --end-idx 150 --init-from-ramp --ramp-len 50 --max-iter 1 --max-wall-sec 300
```

Lemniscate (baseline replay-only subset):
```
python3 examples/run_ilqr_lemniscate.py --backend v1_3 --full-dataset --mode replay_only_full --start-idx 0 --end-idx 150 --max-wall-sec 300
```

Lemniscate (controller with ramp init):
```
python3 examples/run_ilqr_lemniscate.py --backend v1_3 --full-dataset --mode replay_dataset_mode_full --start-idx 0 --end-idx 150 --init-from-ramp --ramp-len 50 --max-iter 1 --max-wall-sec 300
```

## Runtime notes
- Full 320-step iLQR runs can be slow on this machine; use `--start-idx/--end-idx` to run subsets.
- Use `--plot-every` to reduce plot overhead; NPZ is still saved every MPC step.
- If `--max-wall-sec` is hit, partial artifacts and run reports are still written.

## Known limitations
- Full-dataset iLQR may hit max-wall time on slower machines.
- Plotting every step is expensive; use `--plot-every` to throttle.

## GO statement
GO to Phase C.
