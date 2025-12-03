"""
Configuration utilities.
"""

import yaml
from typing import Dict, Any


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to config file
    
    Returns:
        Configuration dictionary
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def save_config(config: Dict[str, Any], config_path: str):
    """
    Save configuration to YAML file.
    
    Args:
        config: Configuration dictionary
        config_path: Path to save config
    """
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)


def get_default_config() -> Dict[str, Any]:
    """
    Get default configuration.
    
    Returns:
        Default config dictionary
    """
    return {
        'catheter': {
            'n_segments': 3,
            'segment_lengths': [0.104, 0.013, 0.029],  # meters
            'material_E': 31.03e6,  # Pa
            'material_G': 8.11e6,   # Pa
        },
        'training': {
            'n_episodes': 5000,
            'batch_size': 256,
            'learning_rate': 3e-4,
            'discount': 0.99,
            'tau': 0.005,
            'exploration_noise': 0.1,
            'start_steps': 10000,
        },
        'environment': {
            'max_steps': 200,
            'dt': 0.05,
            'state_dim': 24,
            'reward_weights': {
                'distance': -1.0,
                'curvature': -0.01,
                'energy': -0.001,
                'smoothness': -0.01,
                'goal_bonus': 10.0,
            }
        }
    }
