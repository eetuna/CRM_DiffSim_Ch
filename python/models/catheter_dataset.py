"""PyTorch dataset for catheter dynamics data."""

import numpy as np
import torch
from torch.utils.data import Dataset
from typing import Optional, Tuple


class CatheterDataset(Dataset):
    """Dataset for catheter state and control data."""
    
    def __init__(
        self,
        states: np.ndarray,
        actions: np.ndarray,
        next_states: np.ndarray,
        normalize: bool = True
    ):
        """Initialize dataset.
        
        Args:
            states: Current states (N, state_dim)
            actions: Control actions (N, action_dim)
            next_states: Next states (N, state_dim)
            normalize: Whether to normalize data
        """
        self.states = torch.FloatTensor(states)
        self.actions = torch.FloatTensor(actions)
        self.next_states = torch.FloatTensor(next_states)
        
        if normalize:
            self.state_mean = self.states.mean(dim=0)
            self.state_std = self.states.std(dim=0) + 1e-8
            self.action_mean = self.actions.mean(dim=0)
            self.action_std = self.actions.std(dim=0) + 1e-8
            
            self.states = (self.states - self.state_mean) / self.state_std
            self.actions = (self.actions - self.action_mean) / self.action_std
            self.next_states = (self.next_states - self.state_mean) / self.state_std
        else:
            self.state_mean = None
            self.state_std = None
            self.action_mean = None
            self.action_std = None
    
    def __len__(self) -> int:
        """Get dataset size."""
        return len(self.states)
    
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Get data sample.
        
        Args:
            idx: Sample index
            
        Returns:
            (state, action, next_state) tuple
        """
        return self.states[idx], self.actions[idx], self.next_states[idx]
    
    def denormalize_state(self, normalized_state: torch.Tensor) -> torch.Tensor:
        """Convert normalized state back to original scale."""
        if self.state_mean is not None:
            return normalized_state * self.state_std + self.state_mean
        return normalized_state
    
    def denormalize_action(self, normalized_action: torch.Tensor) -> torch.Tensor:
        """Convert normalized action back to original scale."""
        if self.action_mean is not None:
            return normalized_action * self.action_std + self.action_mean
        return normalized_action
