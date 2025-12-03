"""
Forward kinematics API - wraps C++ Cosserat BVP/IVP solver.

Based on: "Real-Time Forward Kinematics and Jacobian for Control of an 
MRI-Guided Magnetically Actuated Robotic Catheter" (TMECH 2024)
"""

import numpy as np
from typing import Dict, Optional
import crm_cpp


class ForwardKinematicsAPI:
    """
    Forward kinematics using your C++ BVP/IVP solver.
    
    Wraps CRM::SolveIVP from your implementation.
    """
    
    def __init__(self):
        """Initialize with your CRM solver."""
        self._solver = crm_cpp.ForwardKinematics()
        self._jacobian_computer = crm_cpp.JacobianComputer()
    
    def compute(
        self,
        angles: np.ndarray,
        currents: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Compute forward kinematics.
        
        Args:
            angles: Configuration parameters (N*3 where N=number of segments)
            currents: Coil currents (6D), range [-2, 2] A
        
        Returns:
            Dictionary with:
            - tip_position: (3,) tip position in meters
            - tip_orientation: (9,) rotation matrix flattened
            - positions: (n_points, 3) all centerline positions
            - orientations: (n_points, 9) all orientations
            - success: bool
        """
        angles = np.asarray(angles, dtype=np.float64)
        currents = np.asarray(currents, dtype=np.float64)
        
        assert currents.shape == (6,), f"Expected 6 currents, got {currents.shape}"
        assert np.all(np.abs(currents) <= 2.0), "Currents must be in [-2, 2] A"
        
        # Call your C++ solver
        result = self._solver.compute(angles, currents)
        
        # Extract tip (last point)
        positions = result['positions']
        orientations = result['orientations']
        
        return {
            'tip_position': positions[-1, :],
            'tip_orientation': orientations[-1, :],
            'positions': positions,
            'orientations': orientations,
            'success': result['success']
        }
    
    def get_tip_position(
        self,
        angles: np.ndarray,
        currents: np.ndarray
    ) -> np.ndarray:
        """
        Get only tip position (convenience method).
        
        Args:
            angles: Configuration parameters
            currents: Coil currents (6D)
        
        Returns:
            Tip position (3D) in meters
        """
        angles = np.asarray(angles, dtype=np.float64)
        currents = np.asarray(currents, dtype=np.float64)
        return self._solver.get_tip_position(angles, currents)
    
    def get_tip_orientation(
        self,
        angles: np.ndarray,
        currents: np.ndarray
    ) -> np.ndarray:
        """
        Get tip orientation as rotation matrix.
        
        Args:
            angles: Configuration parameters
            currents: Coil currents (6D)
        
        Returns:
            Rotation matrix (3×3)
        """
        result = self.compute(angles, currents)
        R_flat = result['tip_orientation']
        return R_flat.reshape(3, 3)
    
    def compute_jacobian(
        self,
        angles: np.ndarray,
        currents: np.ndarray
    ) -> np.ndarray:
        """
        Compute Jacobian matrix using your C++ implementation.
        
        Args:
            angles: Configuration parameters
            currents: Coil currents (6D)
        
        Returns:
            Jacobian matrix (6, n_dof)
        """
        angles = np.asarray(angles, dtype=np.float64)
        currents = np.asarray(currents, dtype=np.float64)
        
        assert currents.shape == (6,)
        assert np.all(np.abs(currents) <= 2.0)
        
        return self._jacobian_computer.compute(angles, currents)
    
    def set_params(self, length: float, stiffness: float, damping: float):
        """
        Set catheter parameters.
        
        Args:
            length: Total catheter length (m)
            stiffness: Bending stiffness (N⋅m²)
            damping: Damping coefficient
        """
        self._solver.set_params(length, stiffness, damping)
    
    def batch_compute(
        self,
        angles_batch: np.ndarray,
        currents_batch: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Batch computation for multiple configurations.
        
        Args:
            angles_batch: (batch_size, n_dof)
            currents_batch: (batch_size, 6)
        
        Returns:
            Dictionary with batched results
        """
        batch_size = angles_batch.shape[0]
        
        tip_positions = []
        tip_orientations = []
        successes = []
        
        for i in range(batch_size):
            result = self.compute(angles_batch[i], currents_batch[i])
            tip_positions.append(result['tip_position'])
            tip_orientations.append(result['tip_orientation'])
            successes.append(result['success'])
        
        return {
            'tip_positions': np.array(tip_positions),
            'tip_orientations': np.array(tip_orientations),
            'successes': np.array(successes)
        }


class InverseKinematicsAPI:
    """
    Inverse kinematics using your C++ BVP solver.
    
    Wraps CRM::SolveBVP from your implementation.
    """
    
    def __init__(self):
        """Initialize with your CRM BVP solver."""
        self._solver = crm_cpp.InverseKinematics()
    
    def solve(
        self,
        target_position: np.ndarray,
        currents: np.ndarray,
        target_orientation: Optional[np.ndarray] = None,
        q_init: Optional[np.ndarray] = None
    ) -> Dict[str, np.ndarray]:
        """
        Solve inverse kinematics.
        
        Args:
            target_position: Desired tip position (3D)
            currents: Coil currents (6D)
            target_orientation: Desired orientation (optional, 3×3 or 9D)
            q_init: Initial guess for configuration (optional)
        
        Returns:
            Dictionary with:
            - q: Configuration parameters
            - success: bool
            - iterations: int (BVP solver iterations)
            - residual: float (BVP residual norm)
        """
        target_position = np.asarray(target_position, dtype=np.float64)
        currents = np.asarray(currents, dtype=np.float64)
        
        assert target_position.shape == (3,)
        assert currents.shape == (6,)
        assert np.all(np.abs(currents) <= 2.0)
        
        # Call your C++ BVP solver
        result = self._solver.solve(
            target_position,
            target_orientation,
            currents,
            q_init
        )
        
        return result
    
    def solve_batch(
        self,
        target_positions: np.ndarray,
        currents_batch: np.ndarray
    ) -> Dict[str, np.ndarray]:
        """
        Batch inverse kinematics.
        
        Args:
            target_positions: (batch_size, 3)
            currents_batch: (batch_size, 6)
        
        Returns:
            Dictionary with batched results
        """
        batch_size = target_positions.shape[0]
        
        q_results = []
        successes = []
        
        for i in range(batch_size):
            result = self.solve(target_positions[i], currents_batch[i])
            q_results.append(result['q'])
            successes.append(result['success'])
        
        return {
            'q_batch': np.array(q_results),
            'successes': np.array(successes)
        }
