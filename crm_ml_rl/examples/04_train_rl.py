"""
Train RL agent using hybrid dynamics.
"""

import numpy as np
import torch
import argparse
from pathlib import Path

from crm_ml_rl.api import DynamicsAPI, HybridDynamics
from crm_ml_rl.models import DynamicsModelNN
from crm_ml_rl.rl import MRICatheterEnv, TD3Agent, train_td3
from crm_ml_rl.utils import plot_training_curves


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_episodes', type=int, default=5000)
    parser.add_argument('--state_dim', type=int, default=24)
    parser.add_argument('--dynamics_model', type=str, default=None,
                       help='Path to trained dynamics model (optional)')
    parser.add_argument('--save_path', type=str, default='models/td3_agent.pth')
    args = parser.parse_args()
    
    print("="*60)
    print("Training RL Agent for MRI Catheter Navigation")
    print("="*60)
    
    # Initialize analytical dynamics (your C++ solver)
    print("\nInitializing dynamics...")
    analytical_dynamics = DynamicsAPI(dt=0.05)
    
    # Load ML dynamics model if available
    ml_model = None
    if args.dynamics_model and Path(args.dynamics_model).exists():
        print(f"Loading dynamics model from {args.dynamics_model}")
        ml_model = DynamicsModelNN()
        ml_model.load_state_dict(torch.load(args.dynamics_model))
        ml_model.eval()
        print("  ✓ Dynamics model loaded")
    else:
        print("No dynamics model provided, using analytical only")
    
    # Create hybrid dynamics
    hybrid_dynamics = HybridDynamics(analytical_dynamics, ml_model)
    
    # Create environment
    print("\nCreating environment...")
    env = MRICatheterEnv(
        dynamics=hybrid_dynamics,
        max_steps=200,
        dt=0.05,
        state_dim=args.state_dim
    )
    print(f"  State space: {env.observation_space.shape}")
    print(f"  Action space: {env.action_space.shape}")
    
    # Create agent
    print("\nCreating TD3 agent...")
    agent = TD3Agent(
        state_dim=env.observation_space.shape[0],
        action_dim=env.action_space.shape[0],
        max_action=2.0
    )
    print(f"  Device: {agent.device}")
    
    # Train
    print(f"\n{'='*60}")
    print(f"Starting training for {args.n_episodes} episodes...")
    print(f"{'='*60}\n")
    
    stats = train_td3(
        env=env,
        agent=agent,
        n_episodes=args.n_episodes,
        save_path=args.save_path
    )
    
    # Plot results
    print("\nGenerating training curves...")
    plot_training_curves(stats, save_path='training_curves.png')
    
    print("\n" + "="*60)
    print("Training complete!")
    print(f"Agent saved to {args.save_path}")
    print("="*60)


if __name__ == '__main__':
    main()
