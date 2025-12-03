"""
Tests for dynamics API.
"""

import numpy as np
import pytest
from crm_ml_rl.api import DynamicsAPI


def test_dynamics_init():
    """Test initialization."""
    api = DynamicsAPI(dt=0.05)
    assert api.dt == 0.05
    assert api._simulator is not None


def test_dynamics_step():
    """Test dynamics step."""
    api = DynamicsAPI(dt=0.05)
    
    state = api.reset_state(state_dim=24)
    currents = np.zeros(6)
    
    next_state = api.step(state, currents)
    
    assert next_state.shape == (24,)


def test_dynamics_trajectory():
    """Test trajectory simulation."""
    api = DynamicsAPI(dt=0.05)
    
    initial_state = api.reset_state(state_dim=24)
    n_steps = 10
    control_sequence = np.random.uniform(-1.0, 1.0, size=(n_steps, 6))
    
    trajectory = api.simulate_trajectory(initial_state, control_sequence)
    
    assert trajectory.shape == (n_steps + 1, 24)


def test_timestep():
    """Test timestep get/set."""
    api = DynamicsAPI(dt=0.05)
    
    assert api.get_timestep() == 0.05
    
    api.set_timestep(0.01)
    assert api.get_timestep() == 0.01


if __name__ == '__main__':
    pytest.main([__file__])
