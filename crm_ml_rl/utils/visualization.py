"""
Visualization utilities.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional


def plot_training_curves(
    stats: Dict[str, List[float]],
    save_path: Optional[str] = None
):
    """
    Plot training curves.
    
    Args:
        stats: Training statistics
        save_path: Path to save figure
    """
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # Episode rewards
    axes[0].plot(stats['episode_rewards'], alpha=0.3, label='Episode reward')
    
    # Moving average
    window = 50
    if len(stats['episode_rewards']) >= window:
        moving_avg = np.convolve(
            stats['episode_rewards'],
            np.ones(window) / window,
            mode='valid'
        )
        axes[0].plot(moving_avg, label=f'Moving avg ({window} ep)')
    
    axes[0].set_xlabel('Episode')
    axes[0].set_ylabel('Reward')
    axes[0].set_title('Training Reward')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Evaluation rewards
    if 'eval_rewards' in stats and len(stats['eval_rewards']) > 0:
        eval_episodes = np.linspace(0, len(stats['episode_rewards']),
                                    len(stats['eval_rewards']))
        axes[1].plot(eval_episodes, stats['eval_rewards'], marker='o')
        axes[1].set_xlabel('Episode')
        axes[1].set_ylabel('Reward')
        axes[1].set_title('Evaluation Reward')
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Training curves saved to {save_path}")
    
    plt.show()


def visualize_trajectory(
    positions: np.ndarray,
    target: np.ndarray,
    title: str = "Catheter Trajectory",
    save_path: Optional[str] = None
):
    """
    Visualize 3D catheter trajectory.
    
    Args:
        positions: Trajectory positions (n_steps, 3)
        target: Target position (3,)
        title: Plot title
        save_path: Path to save figure
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot trajectory
    ax.plot(positions[:, 0], positions[:, 1], positions[:, 2],
            'b-', linewidth=2, label='Trajectory')
    
    # Start point
    ax.scatter(positions[0, 0], positions[0, 1], positions[0, 2],
              c='g', s=100, marker='o', label='Start')
    
    # End point
    ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2],
              c='r', s=100, marker='o', label='End')
    
    # Target
    ax.scatter(target[0], target[1], target[2],
              c='gold', s=200, marker='*', label='Target')
    
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Trajectory visualization saved to {save_path}")
    
    plt.show()
