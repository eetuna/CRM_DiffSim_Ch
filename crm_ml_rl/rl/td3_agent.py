"""
TD3 (Twin Delayed Deep Deterministic Policy Gradient) agent.

Based on: "Addressing Function Approximation Error in Actor-Critic Methods"
Fujimoto et al., ICML 2018
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import Tuple
import copy


class Actor(nn.Module):
    """Actor network (policy) for TD3."""
    
    def __init__(
        self,
        state_dim: int,
        action_dim: int,
        hidden_dims: Tuple[int, int] = (256, 256),
        max_action: float = 2.0
    ):
        super().__init__()
        
        self.net = nn.Sequential(
            nn.Linear(state_dim, hidden_dims[0]),
            nn.ReLU(),
            nn.Linear(hidden_dims[0], hidden_dims[1]),
            nn.ReLU(),
            nn.Linear(hidden_dims[1], action_dim),
            nn.Tanh()
        )
        
        self.max_action = max_action
    
    def forward(self, state: torch.Tensor) -> torch.Tensor:
        """Output scaled action in [-max_action, max_action]."""
        return self.max_action * self.net(state)


class Critic(nn.Module):
    """Twin critic networks (Q-functions) for TD3."""
    
    def __init__(
        self,
        state_dim: int,
        action_dim: int,
        hidden_dims: Tuple[int, int] = (256, 256)
    ):
        super().__init__()
        
        # Q1 network
        self.q1 = nn.Sequential(
            nn.Linear(state_dim + action_dim, hidden_dims[0]),
            nn.ReLU(),
            nn.Linear(hidden_dims[0], hidden_dims[1]),
            nn.ReLU(),
            nn.Linear(hidden_dims[1], 1)
        )
        
        # Q2 network
        self.q2 = nn.Sequential(
            nn.Linear(state_dim + action_dim, hidden_dims[0]),
            nn.ReLU(),
            nn.Linear(hidden_dims[0], hidden_dims[1]),
            nn.ReLU(),
            nn.Linear(hidden_dims[1], 1)
        )
    
    def forward(
        self,
        state: torch.Tensor,
        action: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Return both Q-values."""
        sa = torch.cat([state, action], dim=1)
        return self.q1(sa), self.q2(sa)
    
    def Q1(self, state: torch.Tensor, action: torch.Tensor) -> torch.Tensor:
        """Return only Q1 (for actor loss)."""
        sa = torch.cat([state, action], dim=1)
        return self.q1(sa)


class ReplayBuffer:
    """Experience replay buffer for off-policy RL."""
    
    def __init__(
        self,
        state_dim: int,
        action_dim: int,
        max_size: int = 1000000
    ):
        self.max_size = max_size
        self.ptr = 0
        self.size = 0
        
        self.state = np.zeros((max_size, state_dim))
        self.action = np.zeros((max_size, action_dim))
        self.next_state = np.zeros((max_size, state_dim))
        self.reward = np.zeros((max_size, 1))
        self.done = np.zeros((max_size, 1))
    
    def add(
        self,
        state: np.ndarray,
        action: np.ndarray,
        next_state: np.ndarray,
        reward: float,
        done: bool
    ):
        """Add transition to buffer."""
        self.state[self.ptr] = state
        self.action[self.ptr] = action
        self.next_state[self.ptr] = next_state
        self.reward[self.ptr] = reward
        self.done[self.ptr] = done
        
        self.ptr = (self.ptr + 1) % self.max_size
        self.size = min(self.size + 1, self.max_size)
    
    def sample(
        self,
        batch_size: int
    ) -> Tuple[torch.Tensor, ...]:
        """Sample random batch."""
        ind = np.random.randint(0, self.size, size=batch_size)
        
        return (
            torch.FloatTensor(self.state[ind]),
            torch.FloatTensor(self.action[ind]),
            torch.FloatTensor(self.next_state[ind]),
            torch.FloatTensor(self.reward[ind]),
            torch.FloatTensor(self.done[ind])
        )


class TD3Agent:
    """
    TD3 agent for continuous control.
    
    Key features:
    - Twin critics to reduce overestimation
    - Delayed policy updates
    - Target policy smoothing
    """
    
    def __init__(
        self,
        state_dim: int,
        action_dim: int,
        max_action: float = 2.0,
        discount: float = 0.99,
        tau: float = 0.005,
        policy_noise: float = 0.2,
        noise_clip: float = 0.5,
        policy_freq: int = 2,
        lr: float = 3e-4
    ):
        """
        Initialize TD3 agent.
        
        Args:
            state_dim: State dimension
            action_dim: Action dimension
            max_action: Maximum action value
            discount: Discount factor (gamma)
            tau: Target network update rate
            policy_noise: Noise added to target policy
            noise_clip: Clip range for target policy noise
            policy_freq: Frequency of delayed policy updates
            lr: Learning rate
        """
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # Actor
        self.actor = Actor(state_dim, action_dim, max_action=max_action).to(self.device)
        self.actor_target = copy.deepcopy(self.actor)
        self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=lr)
        
        # Critic
        self.critic = Critic(state_dim, action_dim).to(self.device)
        self.critic_target = copy.deepcopy(self.critic)
        self.critic_optimizer = torch.optim.Adam(self.critic.parameters(), lr=lr)
        
        # Hyperparameters
        self.max_action = max_action
        self.discount = discount
        self.tau = tau
        self.policy_noise = policy_noise
        self.noise_clip = noise_clip
        self.policy_freq = policy_freq
        
        self.total_it = 0
    
    def select_action(self, state: np.ndarray) -> np.ndarray:
        """Select action (deterministic, for evaluation)."""
        state = torch.FloatTensor(state.reshape(1, -1)).to(self.device)
        return self.actor(state).cpu().data.numpy().flatten()
    
    def select_action_with_noise(
        self,
        state: np.ndarray,
        noise: float = 0.1
    ) -> np.ndarray:
        """Select action with exploration noise."""
        action = self.select_action(state)
        noise_vec = np.random.normal(0, noise, size=action.shape)
        action = action + noise_vec
        return np.clip(action, -self.max_action, self.max_action)
    
    def train(self, replay_buffer: ReplayBuffer, batch_size: int = 256):
        """Train agent on batch from replay buffer."""
        self.total_it += 1
        
        # Sample batch
        state, action, next_state, reward, done = replay_buffer.sample(batch_size)
        
        state = state.to(self.device)
        action = action.to(self.device)
        next_state = next_state.to(self.device)
        reward = reward.to(self.device)
        done = done.to(self.device)
        
        with torch.no_grad():
            # Target policy smoothing: add clipped noise to next action
            noise = (torch.randn_like(action) * self.policy_noise).clamp(
                -self.noise_clip, self.noise_clip
            )
            
            next_action = (self.actor_target(next_state) + noise).clamp(
                -self.max_action, self.max_action
            )
            
            # Compute target Q-value (minimum of twin Q-functions)
            target_Q1, target_Q2 = self.critic_target(next_state, next_action)
            target_Q = torch.min(target_Q1, target_Q2)
            target_Q = reward + (1 - done) * self.discount * target_Q
        
        # Get current Q estimates
        current_Q1, current_Q2 = self.critic(state, action)
        
        # Critic loss (MSE)
        critic_loss = F.mse_loss(current_Q1, target_Q) + F.mse_loss(current_Q2, target_Q)
        
        # Optimize critic
        self.critic_optimizer.zero_grad()
        critic_loss.backward()
        self.critic_optimizer.step()
        
        # Delayed policy updates
        if self.total_it % self.policy_freq == 0:
            # Actor loss: maximize Q1(s, actor(s))
            actor_loss = -self.critic.Q1(state, self.actor(state)).mean()
            
            # Optimize actor
            self.actor_optimizer.zero_grad()
            actor_loss.backward()
            self.actor_optimizer.step()
            
            # Soft update target networks
            for param, target_param in zip(self.critic.parameters(), self.critic_target.parameters()):
                target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)
            
            for param, target_param in zip(self.actor.parameters(), self.actor_target.parameters()):
                target_param.data.copy_(self.tau * param.data + (1 - self.tau) * target_param.data)
    
    def save(self, filename: str):
        """Save agent."""
        torch.save({
            'actor': self.actor.state_dict(),
            'critic': self.critic.state_dict(),
            'actor_optimizer': self.actor_optimizer.state_dict(),
            'critic_optimizer': self.critic_optimizer.state_dict(),
        }, filename)
    
    def load(self, filename: str):
        """Load agent."""
        checkpoint = torch.load(filename, map_location=self.device)
        self.actor.load_state_dict(checkpoint['actor'])
        self.critic.load_state_dict(checkpoint['critic'])
        self.actor_optimizer.load_state_dict(checkpoint['actor_optimizer'])
        self.critic_optimizer.load_state_dict(checkpoint['critic_optimizer'])
        
        self.actor_target = copy.deepcopy(self.actor)
        self.critic_target = copy.deepcopy(self.critic)
