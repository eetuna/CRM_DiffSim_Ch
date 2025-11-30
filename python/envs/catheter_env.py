"""Gymnasium environment for MRI-actuated catheter control."""

import numpy as np
import gymnasium as gym
from gymnasium import spaces
from typing import Optional, Tuple, Dict, Any

from ..crm_wrapper import CatheterModel


class CatheterEnv(gym.Env):
    """Reinforcement learning environment for catheter navigation."""
    
    metadata = {'render_modes': ['human', 'rgb_array'], 'render_fps': 30}
    
    def __init__(
        self,
        num_segments: int = 20,
        max_field_strength: float = 0.5,
        target_tolerance: float = 0.001,
        max_steps: int = 1000,
        render_mode: Optional[str] = None
    ):
        """Initialize catheter environment.
        
        Args:
            num_segments: Number of catheter segments
            max_field_strength: Maximum magnetic field strength (Tesla)
            target_tolerance: Distance tolerance for reaching target (meters)
            max_steps: Maximum episode steps
            render_mode: Rendering mode ('human' or 'rgb_array')
        """
        super().__init__()
        
        self.num_segments = num_segments
        self.max_field_strength = max_field_strength
        self.target_tolerance = target_tolerance
        self.max_steps = max_steps
        self.render_mode = render_mode
        
        # Initialize catheter model
        self.catheter = CatheterModel(num_segments)
        
        # Action space: 3D magnetic field vector
        self.action_space = spaces.Box(
            low=-max_field_strength,
            high=max_field_strength,
            shape=(3,),
            dtype=np.float32
        )
        
        # Observation space: catheter state + target position
        state_dim = self.catheter.state_dimension
        self.observation_space = spaces.Box(
            low=-np.inf,
            high=np.inf,
            shape=(state_dim + 3,),  # state + target
            dtype=np.float32
        )
        
        # Episode state
        self.state = None
        self.target = None
        self.steps = 0
        
    def reset(
        self,
        seed: Optional[int] = None,
        options: Optional[Dict[str, Any]] = None
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Reset environment to initial state.
        
        Args:
            seed: Random seed
            options: Additional options
            
        Returns:
            Initial observation and info dict
        """
        super().reset(seed=seed)
        
        # Initialize straight catheter
        self.state = np.zeros(self.catheter.state_dimension)
        
        # Random target in workspace
        self.target = self.np_random.uniform(
            low=[0.0, -0.05, -0.05],
            high=[0.1, 0.05, 0.05]
        ).astype(np.float32)
        
        self.steps = 0
        
        obs = self._get_observation()
        info = self._get_info()
        
        return obs, info
    
    def step(
        self,
        action: np.ndarray
    ) -> Tuple[np.ndarray, float, bool, bool, Dict[str, Any]]:
        """Execute one environment step.
        
        Args:
            action: Magnetic field control [Bx, By, Bz]
            
        Returns:
            observation, reward, terminated, truncated, info
        """
        # Clip action to valid range
        action = np.clip(action, -self.max_field_strength, self.max_field_strength)
        
        # Solve forward kinematics
        result = self.catheter.solve_forward_kinematics(
            magnetic_field=action
        )
        
        # Update state
        self.state = result['strains'].flatten()
        tip_position = result['positions'][-1, :]
        
        # Compute reward
        distance_to_target = np.linalg.norm(tip_position - self.target)
        reward = -distance_to_target  # Negative distance as reward
        
        # Check termination
        terminated = distance_to_target < self.target_tolerance
        
        # Check truncation
        self.steps += 1
        truncated = self.steps >= self.max_steps
        
        # Get observation and info
        obs = self._get_observation()
        info = self._get_info()
        info['distance_to_target'] = distance_to_target
        info['tip_position'] = tip_position
        
        return obs, reward, terminated, truncated, info
    
    def _get_observation(self) -> np.ndarray:
        """Get current observation."""
        return np.concatenate([self.state, self.target]).astype(np.float32)
    
    def _get_info(self) -> Dict[str, Any]:
        """Get environment info."""
        return {
            'steps': self.steps,
            'target': self.target.copy()
        }
    
    def render(self):
        """Render environment (placeholder)."""
        if self.render_mode == 'human':
            print(f"Step: {self.steps}, Target: {self.target}")
        elif self.render_mode == 'rgb_array':
            # Return RGB array (implement visualization)
            return np.zeros((480, 640, 3), dtype=np.uint8)
