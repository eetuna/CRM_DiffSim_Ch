"""
Training utilities for RL agent.
"""

import numpy as np
import torch
from typing import Optional
from tqdm import tqdm
from crm_ml_rl.rl.td3_agent import TD3Agent, ReplayBuffer
from crm_ml_rl.rl.environment import MRICatheterEnv


def train_td3(
    env: MRICatheterEnv,
    agent: TD3Agent,
    n_episodes: int = 5000,
    max_steps: int = 200,
    batch_size: int = 256,
    start_steps: int = 10000,
    eval_freq: int = 50,
    save_freq: int = 100,
    save_path: str = "models/td3_agent.pth",
    exploration_noise: float = 0.1
) -> dict:
    """
    Train TD3 agent.
    
    Args:
        env: MRI catheter environment
        agent: TD3 agent
        n_episodes: Number of training episodes
        max_steps: Maximum steps per episode
        batch_size: Batch size for training
        start_steps: Random exploration steps before training
        eval_freq: Evaluation frequency (episodes)
        save_freq: Model saving frequency (episodes)
        save_path: Path to save model
        exploration_noise: Exploration noise std
    
    Returns:
        Dictionary with training statistics
    """
    # Replay buffer
    replay_buffer = ReplayBuffer(
        state_dim=env.observation_space.shape[0],
        action_dim=env.action_space.shape[0]
    )
    
    # Training statistics
    episode_rewards = []
    episode_lengths = []
    eval_rewards = []
    
    total_steps = 0
    
    print(f"Training TD3 agent for {n_episodes} episodes...")
    print(f"State dim: {env.observation_space.shape[0]}, Action dim: {env.action_space.shape[0]}")
    print(f"Using device: {agent.device}")
    
    for episode in tqdm(range(n_episodes), desc="Training"):
        state, _ = env.reset()
        episode_reward = 0
        episode_length = 0
        
        for step in range(max_steps):
            # Select action
            if total_steps < start_steps:
                # Random exploration
                action = env.action_space.sample()
            else:
                # Policy with noise
                action = agent.select_action_with_noise(state, exploration_noise)
            
            # Execute action
            next_state, reward, terminated, truncated, info = env.step(action)
            done = terminated or truncated
            
            # Store transition
            replay_buffer.add(state, action, next_state, reward, done)
            
            state = next_state
            episode_reward += reward
            episode_length += 1
            total_steps += 1
            
            # Train agent
            if total_steps >= start_steps:
                agent.train(replay_buffer, batch_size)
            
            if done:
                break
        
        episode_rewards.append(episode_reward)
        episode_lengths.append(episode_length)
        
        # Evaluation
        if (episode + 1) % eval_freq == 0:
            eval_reward = evaluate_agent(env, agent, n_eval_episodes=10)
            eval_rewards.append(eval_reward)
            
            print(f"\nEpisode {episode + 1}/{n_episodes}")
            print(f"  Train reward: {np.mean(episode_rewards[-eval_freq:]):.2f}")
            print(f"  Eval reward: {eval_reward:.2f}")
            print(f"  Steps: {total_steps}")
        
        # Save model
        if (episode + 1) % save_freq == 0:
            agent.save(save_path)
            print(f"  Model saved to {save_path}")
    
    return {
        'episode_rewards': episode_rewards,
        'episode_lengths': episode_lengths,
        'eval_rewards': eval_rewards
    }


def evaluate_agent(
    env: MRICatheterEnv,
    agent: TD3Agent,
    n_eval_episodes: int = 10
) -> float:
    """
    Evaluate agent performance.
    
    Args:
        env: Environment
        agent: Agent to evaluate
        n_eval_episodes: Number of evaluation episodes
    
    Returns:
        Average episode reward
    """
    eval_rewards = []
    
    for _ in range(n_eval_episodes):
        state, _ = env.reset()
        episode_reward = 0
        done = False
        
        while not done:
            action = agent.select_action(state)  # Deterministic
            state, reward, terminated, truncated, _ = env.step(action)
            episode_reward += reward
            done = terminated or truncated
        
        eval_rewards.append(episode_reward)
    
    return np.mean(eval_rewards)
