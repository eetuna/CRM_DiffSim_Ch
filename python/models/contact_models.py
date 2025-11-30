"""Neural network models for contact force prediction."""

import torch
import torch.nn as nn
from typing import Tuple


class ContactForceNet(nn.Module):
    """Neural network for predicting contact forces."""
    
    def __init__(
        self,
        state_dim: int,
        hidden_dims: Tuple[int, ...] = (256, 256, 128),
        output_dim: int = 3
    ):
        """Initialize contact force network.
        
        Args:
            state_dim: Dimension of catheter state
            hidden_dims: Hidden layer dimensions
            output_dim: Output dimension (force vector)
        """
        super().__init__()
        
        layers = []
        prev_dim = state_dim
        
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.ReLU(),
                nn.LayerNorm(hidden_dim)
            ])
            prev_dim = hidden_dim
        
        layers.append(nn.Linear(prev_dim, output_dim))
        
        self.network = nn.Sequential(*layers)
    
    def forward(self, state: torch.Tensor) -> torch.Tensor:
        """Forward pass.
        
        Args:
            state: Catheter state (batch_size, state_dim)
            
        Returns:
            Predicted contact force (batch_size, output_dim)
        """
        return self.network(state)
