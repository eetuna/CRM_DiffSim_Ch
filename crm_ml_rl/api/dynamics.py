"""
Dynamics API - wraps C++ Hybrid CRDM solver.

Based on: "Free-Space Dynamic Modeling of an MRI-Actuated Robotic Catheter" 
(TMECH 2025)
"""

import numpy as np
from typing import Optional
import crm_cpp


class DynamicsAPI:
    """
    Dynamics simulation using your C++ Hybrid CRDM.
    
    From your dynamics paper (Table IV):
    - Runtime: 6.2ms per step
    - Accuracy: 1-3mm over 10-100ms step sizes
    - Better numerical stability than full-body CRDM
    
    Wraps CRMDYN::IntegrateStep from your implementation.
    """
    
    def __init__(
        self,
        dt: float = 0.05  # 50ms (your 20 Hz servo rate)
    ):
        """
        Initialize dynamics simulator.
        
        Args:
            dt: Time step in seconds (default: 0.05s = 20 Hz)
        """
        self._simulator = crm_cpp.DynamicsSimulator(dt)
        self.dt = dt
    
    def step(
        self,
        state: np.ndarray,
        currents: np.ndarray
    ) -> np.ndarray:
        """
        Step dynamics forward by dt using your CRMDYN solver.
        
        Args:
            state: Current state (variable dimension based on your CRMDYN config)
            currents: Coil currents (6D), range [-2, 2] A
        
        Returns:
            Next state (same dimension as input) after time dt
        """
        state = np.asarray(state, dtype=np.float64)
        currents = np.asarray(currents, dtype=np.float64)
        
        assert currents.shape == (6,), f"Expected 6 currents, got {currents.shape}"
        assert np.all(np.abs(currents) <= 2.0), "Currents must be in [-2, 2] A"
        
        # Call your C++ dynamics integrator
        next_state = self._simulator.step(state, currents)
        
        return next_state
    
    def simulate_trajectory(
        self,
        initial_state: np.ndarray,
        control_sequence: np.ndarray
    ) -> np.ndarray:
        """
        Simulate trajectory using your CRMDYN simulator.
        
        Args:
            initial_state: Starting state (variable dimension)
            control_sequence: Control inputs (n_steps, 6)
        
        Returns:
            State trajectory (n_steps+1, state_dim)
        """
        initial_state = np.asarray(initial_state, dtype=np.float64)
        control_sequence = np.asarray(control_sequence, dtype=np.float64)
        
        assert control_sequence.ndim == 2
        assert control_sequence.shape[1] == 6
        
        # Use your C++ batch simulator
        trajectory = self._simulator.simulate(initial_state, control_sequence)
        
        return trajectory
    
    def reset_state(
        self,
        state_dim: Optional[int] = None,
        initial_position: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Create initial state for simulation.
        
        Args:
            state_dim: State dimension (if known)
            initial_position: Initial tip position (optional)
        
        Returns:
            Initial state vector
        """
        if state_dim is None:
            # Default dimension (adjust based on your CRMDYN configuration)
            state_dim = 24  # Example: [p(3), R(9), u(6), v(3), w(3)]
        
        state = np.zeros(state_dim)
        
        if initial_position is not None:
            state[0:3] = initial_position
        
        return state
    
    def set_timestep(self, dt: float):
        """
        Set time step for dynamics integration.
        
        Args:
            dt: Time step in seconds
        """
        self._simulator.set_timestep(dt)
        self.dt = dt
    
    def get_timestep(self) -> float:
        """
        Get current time step.
        
        Returns:
            Time step in seconds
        """
        return self._simulator.get_timestep()
