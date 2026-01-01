# v1.3 PR-3 Completion Report

## API changes
- Added public entry point `crm_step_v1_3(x, u, *, Li, cfg, return_cache=False)` in `python/crm_diffsims/dynamics/__init__.py`.
- Added v1.3 linearization helper `linearize_v1_3` with operator-first API and optional dense matrices.

## Integration points
- `python/crm_diffsims/dynamics/linearize_v1_3.py`
- `python/crm_diffsims/control/ilqr.py` supports `linearization_backend="v1_3"`.

## Benchmark results
- Autograd Jacobian (full A,B): 77.450506s
- v1.3 operator apply (A·dx, B·du): 0.298129s
- Speedup: 259.8x

Outputs:
- `output_data/v1_3_benchmark.json`
- `output_data/v1_3_benchmark.txt`

## Commands run
- python3 examples/warmup_crm_extension.py
- python3 scripts/benchmark_v1_3_vs_autograd.py
- python3 -m pytest tests/test_v1_3_performance_smoke.py -s

## GO / NO-GO for Phase B controllers
- GO: v1.3 VJP is callable from Python, linearization operators work on stable inputs, and operator-based speedup exceeds 50x vs autograd.
