import torch
import torch.nn as nn
import numpy as np
from collections import deque
import crm_cpp

class ReplayBuffer:
    def __init__(self, capacity: int = 100000):
        self.buffer = deque(maxlen=capacity)
    
    def push(self, state, action, reward, next_state, done):
        self.buffer.append((state, action, reward, next_state, done))
    
    def sample(self, batch_size: int):
        indices = np.random.choice(len(self.buffer), batch_size, replace=False)
        states, actions, rewards, next_states, dones = zip(*[self.buffer[i] for i in indices])
        return (np.array(states), np.array(actions), np.array(rewards),
                np.array(next_states), np.array(dones))
    
    def __len__(self):
        return len(self.buffer)


class Actor(nn.Module):
    """Policy network: state → currents."""
    def __init__(self, state_dim: int = 3, action_dim: int = 3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, action_dim),
            nn.Tanh()  # Output in [-1, 1]
        )
    
    def forward(self, x):
        return self.net(x)


class Critic(nn.Module):
    """Q-network: (state, action) → Q-value."""
    def __init__(self, state_dim: int = 3, action_dim: int = 3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim + action_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 256),
            nn.ReLU(),
            nn.Linear(256, 1)
        )
    
    def forward(self, state, action):
        x = torch.cat([state, action], dim=-1)
        return self.net(x)


class TD3Agent:
    """Twin Delayed DDPG (TD3) for continuous control."""
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3, 
                 device: str = 'cpu', action_scale: float = 0.1):
        self.device = device
        self.action_scale = action_scale
        
        self.actor = Actor(state_dim, action_dim).to(device)
        self.actor_target = Actor(state_dim, action_dim).to(device)
        self.actor_target.load_state_dict(self.actor.state_dict())
        
        self.critic1 = Critic(state_dim, action_dim).to(device)
        self.critic1_target = Critic(state_dim, action_dim).to(device)
        self.critic1_target.load_state_dict(self.critic1.state_dict())
        
        self.critic2 = Critic(state_dim, action_dim).to(device)
        self.critic2_target = Critic(state_dim, action_dim).to(device)
        self.critic2_target.load_state_dict(self.critic2.state_dict())
        
        self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=1e-4)
        self.critic1_optimizer = torch.optim.Adam(self.critic1.parameters(), lr=1e-3)
        self.critic2_optimizer = torch.optim.Adam(self.critic2.parameters(), lr=1e-3)
        
        self.replay_buffer = ReplayBuffer()
        self.update_counter = 0
    
    def select_action(self, state: np.ndarray, noise_scale: float = 0.1) -> np.ndarray:
        """Select action with exploration noise."""
        state_tensor = torch.from_numpy(state).float().unsqueeze(0).to(self.device)
        with torch.no_grad():
            action = self.actor(state_tensor).cpu().numpy()[0]
        noise = np.random.randn(len(action)) * noise_scale
        action = np.clip(action + noise, -1, 1) * self.action_scale
        return action
    
    def train_step(self, batch_size: int = 64, gamma: float = 0.99):
        """One training step."""
        if len(self.replay_buffer) < batch_size:
            return
        
        states, actions, rewards, next_states, dones = self.replay_buffer.sample(batch_size)
        
        states = torch.from_numpy(states).float().to(self.device)
        actions = torch.from_numpy(actions).float().to(self.device)
        rewards = torch.from_numpy(rewards).float().unsqueeze(1).to(self.device)
        next_states = torch.from_numpy(next_states).float().to(self.device)
        dones = torch.from_numpy(dones).float().unsqueeze(1).to(self.device)
        
        # Update critics
        with torch.no_grad():
            next_actions = self.actor_target(next_states)
            noise = torch.randn_like(next_actions) * 0.2
            next_actions = torch.clamp(next_actions + noise, -1, 1)
            
            q1_target = self.critic1_target(next_states, next_actions)
            q2_target = self.critic2_target(next_states, next_actions)
            q_target = torch.min(q1_target, q2_target)
            td_target = rewards + gamma * (1 - dones) * q_target
        
        q1 = self.critic1(states, actions)
        critic1_loss = nn.MSELoss()(q1, td_target)
        self.critic1_optimizer.zero_grad()
        critic1_loss.backward()
        self.critic1_optimizer.step()
        
        q2 = self.critic2(states, actions)
        critic2_loss = nn.MSELoss()(q2, td_target)
        self.critic2_optimizer.zero_grad()
        critic2_loss.backward()
        self.critic2_optimizer.step()
        
        # Update actor (delayed)
        if self.update_counter % 2 == 0:
            actor_actions = self.actor(states)
            actor_loss = -self.critic1(states, actor_actions).mean()
            self.actor_optimizer.zero_grad()
            actor_loss.backward()
            self.actor_optimizer.step()
        
        # Soft update targets
        tau = 0.005
        for param, target_param in zip(self.critic1.parameters(), 
                                       self.critic1_target.parameters()):
            target_param.data.copy_(tau * param.data + (1 - tau) * target_param.data)
        
        self.update_counter += 1


def train_rl_policy(crm_dyn: crm_cpp.CosseratRodDynamics, num_episodes: int = 100):
    """Train RL policy using C++ simulator."""
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    agent = TD3Agent(state_dim=3, action_dim=3, device=device, action_scale=0.1)
    
    print("=== RL Policy Training ===")
    
    for episode in range(num_episodes):
        crm_dyn.reset()
        state = np.array(crm_dyn.get_state()[:3])  # Tip position
        episode_reward = 0.0
        
        target = np.random.randn(3) * 0.02  # Random target in workspace
        
        for step in range(100):
            action = agent.select_action(state - target, noise_scale=0.05)
            crm_dyn.step(action, dt=0.01)
            next_state = np.array(crm_dyn.get_state()[:3])
            
            # Reward: negative distance to target
            reward = -np.linalg.norm(next_state - target)
            done = (np.linalg.norm(next_state - target) < 1e-3)
            
            agent.replay_buffer.push(state - target, action, reward, 
                                     next_state - target, done)
            agent.train_step()
            
            state = next_state
            episode_reward += reward
            
            if done:
                break
        
        if (episode + 1) % 10 == 0:
            print(f"Episode {episode+1}: Return = {episode_reward:.3f}")
    
    torch.save(agent.actor.state_dict(), "rl_policy.pth")
    print("RL policy saved to rl_policy.pth")


if __name__ == "__main__":
    # Pseudo-code: replace with actual initialization
    # crm_dyn = crm_cpp.CosseratRodDynamics("config.json")
    # train_rl_policy(crm_dyn)
    print("RL training script ready. Uncomment main to run.")
