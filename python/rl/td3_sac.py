"""TD3 and SAC implementations for catheter control."""

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from typing import Tuple, Optional


class Actor(nn.Module):
    """Actor network for continuous control."""
    
    def __init__(self, state_dim: int, action_dim: int, max_action: float):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(state_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, action_dim),
            nn.Tanh()
        )
        self.max_action = max_action
    
    def forward(self, state: torch.Tensor) -> torch.Tensor:
        return self.max_action * self.network(state)


class Critic(nn.Module):
    """Critic network for Q-value estimation."""
    
    def __init__(self, state_dim: int, action_dim: int):
        super().__init__()
        self.q1 = nn.Sequential(
            nn.Linear(state_dim + action_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, 1)
        )
        self.q2 = nn.Sequential(
            nn.Linear(state_dim + action_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, 1)
        )
    
    def forward(
        self,
        state: torch.Tensor,
        action: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        sa = torch.cat([state, action], dim=1)
        return self.q1(sa), self.q2(sa)


class TD3Agent:
    """Twin Delayed DDPG agent."""
    
    def __init__(
        self,
        state_dim: int,
        action_dim: int,
        max_action: float,
        lr: float = 3e-4,
        gamma: float = 0.99,
        tau: float = 0.005
    ):
        self.actor = Actor(state_dim, action_dim, max_action)
        self.actor_target = Actor(state_dim, action_dim, max_action)
        self.actor_target.load_state_dict(self.actor.state_dict())
        self.actor_optimizer = optim.Adam(self.actor.parameters(), lr=lr)
        
        self.critic = Critic(state_dim, action_dim)
        self.critic_target = Critic(state_dim, action_dim)
        self.critic_target.load_state_dict(self.critic.state_dict())
        self.critic_optimizer = optim.Adam(self.critic.parameters(), lr=lr)
        
        self.max_action = max_action
        self.gamma = gamma
        self.tau = tau
        
    def select_action(self, state: np.ndarray) -> np.ndarray:
        state = torch.FloatTensor(state).unsqueeze(0)
        return self.actor(state).detach().cpu().numpy()[0]
    
    def train(self, replay_buffer, batch_size: int = 256):
        # Sample batch
        state, action, next_state, reward, done = replay_buffer.sample(batch_size)
        
        # Critic update
        with torch.no_grad():
            next_action = self.actor_target(next_state)
            target_q1, target_q2 = self.critic_target(next_state, next_action)
            target_q = torch.min(target_q1, target_q2)
            target_q = reward + (1 - done) * self.gamma * target_q
        
        current_q1, current_q2 = self.critic(state, action)
        critic_loss = nn.MSELoss()(current_q1, target_q) + nn.MSELoss()(current_q2, target_q)
        
        self.critic_optimizer.zero_grad()
        critic_loss.backward()
        self.critic_optimizer.step()
        
        # Actor update (delayed)
        actor_loss = -self.critic.q1(state, self.actor(state)).mean()
        
        self.actor_optimizer.zero_grad()
        actor_loss.backward()
        self.actor_optimizer.step()
        
        # Target update
        for param, target_param in zip(self.critic.parameters(), self.critic_target.parameters()):
            target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)
        
        for param, target_param in zip(self.actor.parameters(), self.actor_target.parameters()):
            target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)


class SACAgent:
    """Soft Actor-Critic agent (placeholder structure)."""
    
    def __init__(self, state_dim: int, action_dim: int, max_action: float):
        # SAC implementation
        pass
