# Autodiff Audit Closeout

## Milestone Status
- v1.0: GO
- v1.1: GO (FD bridge)
- v1.2: GO
- Stress: GO

## References
- v1.2 summary: docs/autodiff_audits/v1_2_final_completion_summary.md
- CP-15: docs/autodiff_audits/CP-15_completion.md
- CP-18: docs/autodiff_audits/CP-18_completion.md

## How To Reproduce
### CP-15
1) python3 -m pytest tests/test_crm_step_fd.py -s
2) python3 -m pytest tests/test_crm_step_rollout_smoke.py -s
3) python3 -m pytest tests/test_crm_step_gradcheck.py -s
4) python3 -m pytest tests/test_crm_step_vjp_hybrid.py -s
5) python3 -m pytest tests/test_crm_step_A_step_analytic.py -s
6) python3 -m pytest tests/test_crm_step_B_step.py -s
7) python3 -m pytest tests/test_crm_step_C_step.py -s

### CP-18
1) python3 -m pytest tests/test_stress_replay_npz.py -s
2) python3 -m pytest tests/test_workspace_branching_943.py -s
3) python3 examples/run_workspace_stress_report.py

## Output Locations
- output_data/stress_report_943.csv
- output_data/stress_outliers_943.csv
