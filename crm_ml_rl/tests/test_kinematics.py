"""
Tests for kinematics API.
"""

import numpy as np
import pytest
from crm_ml_rl.api import ForwardKinematicsAPI, InverseKinematicsAPI


def test_forward_kinematics_init():
    """Test initialization."""
    api = ForwardKinematicsAPI()
    assert api._solver is not None
    assert api._jacobian_computer is not None


def test_forward_kinematics_compute():
    """Test forward kinematics computation."""
    api = ForwardKinematicsAPI()
    
    angles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    currents = np.zeros(6)
    
    result = api.compute(angles, currents)
    
    assert 'tip_position' in result
    assert 'tip_orientation' in result
    assert 'success' in result
    assert result['success']
    assert result['tip_position'].shape == (3,)


def test_forward_kinematics_tip_position():
    """Test get tip position."""
    api = ForwardKinematicsAPI()
    
    angles = np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.0])
    currents = np.ones(6) * 0.5
    
    tip_pos = api.get_tip_position(angles, currents)
    assert tip_pos.shape == (3,)


def test_forward_kinematics_jacobian():
    """Test Jacobian computation."""
    api = ForwardKinematicsAPI()
    
    angles = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    currents = np.zeros(6)
    
    J = api.compute_jacobian(angles, currents)
    assert J.shape == (6, len(angles))


def test_inverse_kinematics():
    """Test inverse kinematics."""
    api = InverseKinematicsAPI()
    
    target_pos = np.array([0.05, 0.02, 0.08])
    currents = np.zeros(6)
    
    result = api.solve(target_pos, currents)
    
    assert 'q' in result
    assert 'success' in result


def test_batch_compute():
    """Test batch computation."""
    api = ForwardKinematicsAPI()
    
    batch_size = 5
    angles_batch = np.random.uniform(-0.1, 0.1, size=(batch_size, 6))
    currents_batch = np.random.uniform(-1.0, 1.0, size=(batch_size, 6))
    
    results = api.batch_compute(angles_batch, currents_batch)
    
    assert results['tip_positions'].shape == (batch_size, 3)


if __name__ == '__main__':
    pytest.main([__file__])
