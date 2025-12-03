"""
Hybrid predictors: Analytical (C++) + ML corrections.

Implements:
    p_final = p_analytical + Δp_ML
    state_next = state_analytical + Δstate_ML
"""

import numpy as np
import torch
from typing import Optional
from crm_ml_rl.api.kinematics import ForwardKinematicsAPI
from crm_ml_rl.api.dynamics import DynamicsAPI


class HybridKinematics:
    """
    Hybrid kinematics: Your C++ BVP/IVP solver + ML residual correction.
    
    Prediction = Analytical + ML_correction
    """
    
    def __init__(
        self,
        analytical_solver: ForwardKinematicsAPI,
        ml_model: Optional[torch.nn.Module] = None
    ):
        """
        Initialize hybrid kinematics predictor.
        
        Args:
            analytical_solver: Your C++ forward kinematics
            ml_model: Trained residual correction model (optional)
        """
        self.analytical = analytical_solver
        self.ml_model = ml_model
        
        if ml_model is not None:
            self.ml_model.eval()
    
    def predict(
        self,
        angles: np.ndarray,
        currents: np.ndarray
    ) -> dict:
        """
        Hybrid prediction: analytical + ML correction.
        
        Args:
            angles: Configuration parameters
            currents: Coil currents (6D)
        
        Returns:
            Dictionary with:
            - tip_position_analytical: From C++ solver
            - tip_position_corrected: Analytical + ML correction
            - residual: ML correction (3D)
            - tip_orientation: From C++ solver
        """
        # Step 1: Analytical prediction (your C++ solver)
        result = self.analytical.compute(angles, currents)
        p_analytical = result['tip_position']
        
        # Step 2: ML correction (if model available)
        if self.ml_model is not None:
            with torch.no_grad():
                # Prepare input: [currents, angles]
                x = torch.FloatTensor(np.concatenate([currents, angles]))
                x = x.unsqueeze(0)  # Add batch dimension
                
                # Predict residual
                residual = self.ml_model(x).numpy().squeeze()
        else:
            residual = np.zeros(3)
        
        # Step 3: Combine
        p_corrected = p_analytical + residual
        
        return {
            'tip_position_analytical': p_analytical,
            'tip_position_corrected': p_corrected,
            'residual': residual,
            'tip_orientation': result['tip_orientation']
        }
    
    def predict_batch(
        self,
        angles_batch: np.ndarray,
        currents_batch: np.ndarray
    ) -> dict:
        """
        Batch hybrid prediction.
        
        Args:
            angles_batch: (batch_size, n_dof)
            currents_batch: (batch_size, 6)
        
        Returns:
            Dictionary with batched predictions
        """
        batch_size = angles_batch.shape[0]
        
        # Analytical predictions
        analytical_results = self.analytical.batch_compute(
            angles_batch, currents_batch
        )
        p_analytical = analytical_results['tip_positions']
        
        # ML corrections
        if self.ml_model is not None:
            with torch.no_grad():
                x = torch.FloatTensor(
                    np.concatenate([currents_batch, angles_batch], axis=1)
                )
                residuals = self.ml_model(x).numpy()
        else:
            residuals = np.zeros((batch_size, 3))
        
        # Combine
        p_corrected = p_analytical + residuals
        
        return {
            'tip_positions_analytical': p_analytical,
            'tip_positions_corrected': p_corrected,
            'residuals': residuals
        }


class HybridDynamics:
    """
    Hybrid dynamics: Your C++ Hybrid CRDM + ML correction.
    
    State_next = Analytical + ML_correction
    """
    
    def __init__(
        self,
        analytical_simulator: DynamicsAPI,
        ml_model: Optional[torch.nn.Module] = None
    ):
        """
        Initialize hybrid dynamics.
        
        Args:
            analytical_simulator: Your C++ hybrid CRDM
            ml_model: Trained dynamics correction model (optional)
        """
        self.analytical = analytical_simulator
        self.ml_model = ml_model
        
        if ml_model is not None:
            self.ml_model.eval()
    
    def step(
        self,
        state: np.ndarray,
        currents: np.ndarray
    ) -> dict:
        """
        Hybrid dynamics step.
        
        Args:
            state: Current state
            currents: Control input (6D)
        
        Returns:
            Dictionary with analytical and corrected predictions
        """
        # Step 1: Analytical dynamics (your C++ hybrid CRDM)
        state_analytical = self.analytical.step(state, currents)
        
        # Step 2: ML correction (if model available)
        if self.ml_model is not None:
            with torch.no_grad():
                # Simplified state for ML
                state_simplified = self._simplify_state(state)
                x = torch.FloatTensor(
                    np.concatenate([state_simplified, currents])
                )
                x = x.unsqueeze(0)
                
                # Predict correction
                correction_simplified = self.ml_model(x).numpy().squeeze()
                
                # Map back to full state
                correction = self._expand_correction(
                    correction_simplified,
                    state.shape[0]
                )
        else:
            correction = np.zeros_like(state)
        
        # Step 3: Combine
        state_corrected = state_analytical + correction
        
        return {
            'state_analytical': state_analytical,
            'state_corrected': state_corrected,
            'correction': correction
        }
    
    def _simplify_state(self, state: np.ndarray) -> np.ndarray:
        """
        Simplify full state to ML-friendly representation.
        
        This is a placeholder - adjust based on your state structure.
        Example: Extract key components for ML (position, velocities, etc.)
        """
        # Default: use first 18 components
        return state[:18] if len(state) >= 18 else state
    
    def _expand_correction(
        self,
        correction_simplified: np.ndarray,
        full_dim: int
    ) -> np.ndarray:
        """
        Expand ML correction to full state dimension.
        
        This is a placeholder - adjust based on your state structure.
        """
        correction_full = np.zeros(full_dim)
        
        # Map simplified correction to full state
        n = min(len(correction_simplified), full_dim)
        correction_full[:n] = correction_simplified[:n]
        
        return correction_full
