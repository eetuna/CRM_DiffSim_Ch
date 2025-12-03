"""
Tests for RL components.
"""

import numpy as np
import pytest
from crm_ml_rl.api import DynamicsAPI, HybridDynamics
from crm_ml_rl.rl import MRICatheterEnv, TD3Agent


def test_environment_init():
    """Test environment initialization."""
    dynamics = DynamicsAPI(dt=0.05)
    hybrid = HybridDynamics(dynamics, ml_model=None)
    
    env = MRICatheterEnv(dynamics=hybrid, state_dim=24)
    
    assert env.observation_space.shape == (18,)
    assert env.action_space.shape == (6,)


def test_environment_reset():
    """Test environment reset."""
    dynamics = DynamicsAPI(dt=0.05)
    hybrid = HybridDynamics(dynamics, ml_model=None)
    env = MRICatheterEnv(dynamics=hybrid, state_dim=24)
    
    state, info = env.reset()
    
    assert state.shape == (18,)
    assert 'target' in info
    assert 'distance_to_goal' in info


def test_environment_step():
    """Test environment step."""
    dynamics = DynamicsAPI(dt=0.05)
    hybrid = HybridDynamics(dynamics, ml_model=None)
    env = MRICatheterEnv(dynamics=hybrid, state_dim=24)
    
    state, _ = env.reset()
    action = env.action_space.sample()
    
    next_state, reward, terminated, truncated, info = env.step(action)
    
    assert next_state.shape == (18,)
    assert isinstance(reward, (float, np.floating))
    assert isinstance(terminated, (bool, np.bool_))
    assert isinstance(truncated, (bool, np.bool_))
    assert 'distance_to_goal' in info


def test_td3_agent_init():
    """Test TD3 agent initialization."""
    agent = TD3Agent(state_dim=18, action_dim=6, max_action=2.0)
    
    assert agent.max_action == 2.0
    assert agent.actor is not None
    assert agent.critic is not None


def test_td3_agent_action():
    """Test action selection."""
    agent = TD3Agent(state_dim=18, action_dim=6, max_action=2.0)
    
    state = np.random.randn(18)
    action = agent.select_action(state)
    
    assert action.shape == (6,)
    assert np.all(np.abs(action) <= 2.0)


def test_td3_agent_save_load(tmp_path):
    """Test agent save/load."""
    agent = TD3Agent(state_dim=18, action_dim=6, max_action=2.0)
    
    save_path = tmp_path / "agent.pth"
    agent.save(str(save_path))
    
    agent2 = TD3Agent(state_dim=18, action_dim=6, max_action=2.0)
    agent2.load(str(save_path))
    
    # Check that loaded agent produces same output
    state = np.random.randn(18)
    action1 = agent.select_action(state)
    action2 = agent2.select_action(state)
    
    np.testing.assert_array_almost_equal(action1, action2)


if __name__ == '__main__':
    pytest.main([__file__])
