"""
Gymnasium environment for MRI-actuated catheter control.

Uses your hybrid CRDM for dynamics simulation.
"""

import gymnasium as gym
import numpy as np
from typing import Optional, Tuple, Dict, Any
from crm_ml_rl.api.hybrid import HybridDynamics


class MRICatheterEnv(gym.Env):
    """
    Gymnasium environment for MRI catheter navigation.
    
    State space (18D simplified):
        - Tip position: p ∈ ℝ³
        - Tip orientation: R ∈ ℝ³ (simplified)
        - Curvatures: κ ∈ ℝ⁶
        - Linear velocity: v ∈ ℝ³
        - Angular velocity: ω ∈ ℝ³
    
    Action space (6D):
        - Coil currents: I ∈ [-2, 2]⁶ Amps
    
    Uses your hybrid CRDM (C++ + ML) for physics simulation.
    """
    
    metadata = {'render_modes': ['human', 'rgb_array'], 'render_fps': 20}
    
    def __init__(
        self,
        dynamics: HybridDynamics,
        target_position: Optional[np.ndarray] = None,
        max_steps: int = 200,
        dt: float = 0.05,
        reward_weights: Optional[Dict[str, float]] = None,
        state_dim: int = 24  # Your CRMDYN state dimension
    ):
        """
        Initialize environment.
        
        Args:
            dynamics: Hybrid dynamics model (analytical + ML)
            target_position: Goal position (3D), default: random
            max_steps: Maximum episode steps
            dt: Time step (default: 0.05s = 20 Hz)
            reward_weights: Reward function weights
            state_dim: Full state dimension from your CRMDYN
        """
        super().__init__()
        
        self.dynamics = dynamics
        self.max_steps = max_steps
        self.dt = dt
        self.state_dim = state_dim
        
        # Observation space (18D simplified for RL)
        self.observation_space = gym.spaces.Box(
            low=-np.inf,
            high=np.inf,
            shape=(18,),
            dtype=np.float32
        )
        
        # Action space (6D coil currents)
        self.action_space = gym.spaces.Box(
            low=-2.0,
            high=2.0,
            shape=(6,),
            dtype=np.float32
        )
        
        # Reward weights
        if reward_weights is None:
            reward_weights = {
                'distance': -1.0,
                'curvature': -0.01,
                'energy': -0.001,
                'smoothness': -0.01,
                'goal_bonus': 10.0
            }
        self.reward_weights = reward_weights
        
        # Safety constraints (from your prototype)
        self.I_max = 2.0  # Max current (A)
        self.kappa_max = 0.05  # Max curvature (rad/mm)
        
        # Internal state
        self.state_full = None  # Full state (state_dim)
        self.state_obs = None  # Observation state (18D)
        self.target = target_position
        self.step_count = 0
        self.prev_action = np.zeros(6)
    
    def reset(
        self,
        seed: Optional[int] = None,
        options: Optional[dict] = None
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Reset environment."""
        super().reset(seed=seed)
        
        # Reset state to initial configuration
        self.state_full = np.zeros(self.state_dim)
        
        # Initialize orientation (identity matrix flattened, if present in state)
        if self.state_dim >= 12:
            # Assuming state structure: [p(3), R(9), ...]
            self.state_full[3:12] = np.eye(3).flatten()
        
        self.state_obs = self._full_to_obs_state(self.state_full)
        
        # Set random target if not specified
        if self.target is None:
            self.target = self._sample_target()
        
        self.step_count = 0
        self.prev_action = np.zeros(6)
        
        info = {
            'target': self.target,
            'distance_to_goal': self._compute_distance_to_goal()
        }
        
        return self.state_obs.astype(np.float32), info
    
    def step(
        self,
        action: np.ndarray
    ) -> Tuple[np.ndarray, float, bool, bool, Dict[str, Any]]:
        """Execute one step."""
        # Clip action to valid range
        action = np.clip(action, -self.I_max, self.I_max)
        
        # Step dynamics (uses your hybrid CRDM)
        result = self.dynamics.step(self.state_full, action)
        self.state_full = result['state_corrected']
        self.state_obs = self._full_to_obs_state(self.state_full)
        
        # Compute reward
        reward = self._compute_reward(action)
        
        # Check termination
        self.step_count += 1
        distance = self._compute_distance_to_goal()
        
        terminated = distance < 0.005  # 5mm goal tolerance
        truncated = self.step_count >= self.max_steps
        
        # Info
        info = {
            'distance_to_goal': distance,
            'curvature': self._get_curvature_magnitude(),
            'energy': np.sum(action**2),
            'step': self.step_count
        }
        
        self.prev_action = action
        
        return self.state_obs.astype(np.float32), reward, terminated, truncated, info
    
    def _compute_reward(self, action: np.ndarray) -> float:
        """
        Compute reward based on task progress and physics.
        
        Reward = distance_term + curvature_penalty + energy_penalty + 
                 smoothness_penalty + goal_bonus
        """
        w = self.reward_weights
        
        # 1. Distance to goal (primary objective)
        distance = self._compute_distance_to_goal()
        r_distance = w['distance'] * distance
        
        # 2. Curvature penalty (safety constraint)
        curvature_mag = self._get_curvature_magnitude()
        curvature_penalty = w['curvature'] * curvature_mag**2
        
        # 3. Energy penalty (minimize magnetic work)
        energy = np.sum(action**2)
        r_energy = w['energy'] * energy
        
        # 4. Smoothness penalty (avoid jerky motions)
        action_change = np.linalg.norm(action - self.prev_action)
        r_smoothness = w['smoothness'] * action_change
        
        # 5. Goal bonus
        goal_bonus = w['goal_bonus'] if distance < 0.005 else 0.0
        
        total_reward = (r_distance + curvature_penalty + 
                       r_energy + r_smoothness + goal_bonus)
        
        return total_reward
    
    def _compute_distance_to_goal(self) -> float:
        """Compute Euclidean distance from tip to target."""
        tip_position = self.state_obs[0:3]
        return np.linalg.norm(tip_position - self.target)
    
    def _get_curvature_magnitude(self) -> float:
        """Get magnitude of curvatures."""
        # Assuming curvatures are at indices 6:12 in observation
        if len(self.state_obs) >= 12:
            kappa = self.state_obs[6:12]
            return np.linalg.norm(kappa)
        return 0.0
    
    def _full_to_obs_state(self, state_full: np.ndarray) -> np.ndarray:
        """
        Convert full state to observation state (18D).
        
        This is a simplified mapping. Adjust based on your actual state structure.
        
        Assumed full state structure (example):
        [p(3), R(9), u(6), v(3), w(3), ...] = 24+D
        
        Simplified obs state:
        [p(3), R_simplified(3), u(6), v(3), w(3)] = 18D
        """
        obs = np.zeros(18)
        
        # Position (first 3 components)
        obs[0:3] = state_full[0:3]
        
        # Orientation (simplified to first column of R if present)
        if len(state_full) >= 12:
            R_flat = state_full[3:12]
            R = R_flat.reshape(3, 3)
            obs[3:6] = R[:, 0]  # First column
        
        # Curvatures (if present)
        if len(state_full) >= 18:
            obs[6:12] = state_full[12:18]
        
        # Velocities (if present)
        if len(state_full) >= 24:
            obs[12:15] = state_full[18:21]  # Linear velocity
            obs[15:18] = state_full[21:24]  # Angular velocity
        
        return obs
    
    def _sample_target(self) -> np.ndarray:
        """Sample random reachable target position."""
        # Sample within workspace (approximate)
        # Your catheter length: 146mm = 0.146m
        radius = 0.10  # 100mm reach
        
        theta = self.np_random.uniform(0, 2 * np.pi)
        r = self.np_random.uniform(0.03, radius)
        z = self.np_random.uniform(-0.05, 0.05)
        
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        
        return np.array([x, y, z])
    
    def render(self):
        """Render environment (placeholder)."""
        pass
    
    def close(self):
        """Clean up resources."""
        pass
