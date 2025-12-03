"""
Evaluate trained RL agent.
"""

import numpy as np
import torch
import argparse
from pathlib import Path

from crm_ml_rl.api import DynamicsAPI, HybridDynamics
from crm_ml_rl.models import DynamicsModelNN
from crm_ml_rl.rl import MRICatheterEnv, TD3Agent
from crm_ml_rl.utils import visualize_trajectory


def evaluate_agent(
    env: MRICatheterEnv,
    agent: TD3Agent,
    n_trials: int = 100
):
    """Evaluate agent performance."""
    
    successes = 0
    distances = []
    episode_lengths = []
    
    print(f"Evaluating agent for {n_trials} trials...")
    
    for trial in range(n_trials):
        state, info = env.reset()
        positions = [state[0:3]]
        done = False
        steps = 0
        
        while not done:
            action = agent.select_action(state)
            state, reward, terminated, truncated, info = env.step(action)
            positions.append(state[0:3])
            done = terminated or truncated
            steps += 1
        
        # Check success
        final_distance = info['distance_to_goal']
        distances.append(final_distance)
        episode_lengths.append(steps)
        
        if final_distance < 0.005:  # 5mm threshold
            successes += 1
        
        # Visualize first trial
        if trial == 0:
            positions = np.array(positions)
            visualize_trajectory(
                positions,
                env.target,
                title=f"Trial 1: Distance = {final_distance*1000:.1f}mm",
                save_path='evaluation_trajectory.png'
            )
    
    # Print statistics
    print("\n" + "="*60)
    print("EVALUATION RESULTS")
    print("="*60)
    print(f"Success rate: {successes}/{n_trials} ({100*successes/n_trials:.1f}%)")
    print(f"Mean distance: {np.mean(distances)*1000:.2f} ± {np.std(distances)*1000:.2f} mm")
    print(f"Median distance: {np.median(distances)*1000:.2f} mm")
    print(f"Min distance: {np.min(distances)*1000:.2f} mm")
    print(f"Max distance: {np.max(distances)*1000:.2f} mm")
    print(f"Mean episode length: {np.mean(episode_lengths):.1f} steps")
    print("="*60)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--agent', type=str, required=True,
                       help='Path to trained agent')
    parser.add_argument('--dynamics_model', type=str, default=None)
    parser.add_argument('--n_trials', type=int, default=100)
    parser.add_argument('--state_dim', type=int, default=24)
    args = parser.parse_args()
    
    print("="*60)
    print("Evaluating RL Agent")
    print("="*60)
    
    # Initialize environment
    print("\nInitializing environment...")
    analytical_dynamics = DynamicsAPI(dt=0.05)
    
    ml_model = None
    if args.dynamics_model and Path(args.dynamics_model).exists():
        print(f"Loading dynamics model from {args.dynamics_model}")
        ml_model = DynamicsModelNN()
        ml_model.load_state_dict(torch.load(args.dynamics_model))
        ml_model.eval()
    
    hybrid_dynamics = HybridDynamics(analytical_dynamics, ml_model)
    
    env = MRICatheterEnv(
        dynamics=hybrid_dynamics,
        max_steps=200,
        dt=0.05,
        state_dim=args.state_dim
    )
    
    # Load agent
    print(f"\nLoading agent from {args.agent}...")
    agent = TD3Agent(
        state_dim=env.observation_space.shape[0],
        action_dim=env.action_space.shape[0],
        max_action=2.0
    )
    agent.load(args.agent)
    print("  ✓ Agent loaded")
    
    # Evaluate
    print()
    evaluate_agent(env, agent, args.n_trials)


if __name__ == '__main__':
    main()
