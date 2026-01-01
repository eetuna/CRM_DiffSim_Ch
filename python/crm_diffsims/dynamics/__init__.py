"""Dynamics API surface."""

from crm_diffsims.dynamics.step import crm_step
from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3 as _crm_step_v1_3


def crm_step_v1_3(x_t, u_t, *, Li, cfg, return_cache=False):
    if return_cache:
        from crm_diffsims.dynamics.step_v1_3 import crm_step_v1_3_forward

        return crm_step_v1_3_forward(x_t, u_t, Li, cfg)
    return _crm_step_v1_3(x_t, u_t, Li, cfg)

__all__ = ["crm_step", "crm_step_v1_3"]
