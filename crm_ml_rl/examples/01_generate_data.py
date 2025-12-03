"""
Generate training data for kinematics and dynamics models.

This script uses your C++ analytical solvers to generate data.
NOTE: For real residual learning, you need experimental measurements!
"""

import numpy as np
import argparse
from pathlib import Path
from tqdm import tqdm
from crm_ml_rl.api import ForwardKinematicsAPI, DynamicsAPI


def generate_kinematics_data(
    n_samples: int,
    output_path: str,
    config_dim: int = 6  # Adjust based on your catheter config
):
    """
    Generate kinematics training data.
    
    For residual learning, you would need REAL measurements.
    This generates analytical data only (for testing pipeline).
    """
    print(f"Generating {n_samples} kinematics samples...")
    
    # Initialize API (your C++ solver)
    fk_api = ForwardKinematicsAPI()
    
    # Storage
    angles_data = []
    currents_data = []
    positions_data = []
    orientations_data = []
    
    for _ in tqdm(range(n_samples)):
        # Sample random configurations
        angles = np.random.uniform(-np.pi/4, np.pi/4, size=config_dim)
        currents = np.random.uniform(-2.0, 2.0, size=6)
        
        # Compute kinematics
        try:
            result = fk_api.compute(angles, currents)
            if result['success']:
                angles_data.append(angles)
                currents_data.append(currents)
                positions_data.append(result['tip_position'])
                orientations_data.append(result['tip_orientation'])
        except Exception as e:
            print(f"Error: {e}")
            continue
    
    # Save
    np.savez(
        output_path,
        angles=np.array(angles_data),
        currents=np.array(currents_data),
        positions=np.array(positions_data),
        orientations=np.array(orientations_data)
    )
    
    print(f"Saved {len(angles_data)} samples to {output_path}")
    print(f"NOTE: For residual learning, replace with experimental data!")


def generate_dynamics_data(
    n_trajectories: int,
    trajectory_length: int,
    output_path: str,
    state_dim: int = 24  # Adjust based on your CRMDYN
):
    """
    Generate dynamics training data.
    
    For residual learning, you would need REAL measurements.
    This generates analytical data only (for testing pipeline).
    """
    print(f"Generating {n_trajectories} dynamics trajectories...")
    
    # Initialize API (your C++ simulator)
    dynamics_api = DynamicsAPI(dt=0.05)
    
    # Storage
    states_data = []
    actions_data = []
    next_states_data = []
    
    for _ in tqdm(range(n_trajectories)):
        # Initial state
        state = dynamics_api.reset_state(state_dim=state_dim)
        
        for _ in range(trajectory_length):
            # Random action
            action = np.random.uniform(-2.0, 2.0, size=6)
            
            # Step dynamics
            try:
                next_state = dynamics_api.step(state, action)
                
                # Store (use simplified state if needed)
                states_data.append(state[:18])  # Simplified
                actions_data.append(action)
                next_states_data.append(next_state[:18])  # Simplified
                
                state = next_state
            except Exception as e:
                print(f"Error: {e}")
                break
    
    # Save
    np.savez(
        output_path,
        states=np.array(states_data),
        actions=np.array(actions_data),
        next_states=np.array(next_states_data)
    )
    
    print(f"Saved {len(states_data)} samples to {output_path}")
    print(f"NOTE: For residual learning, replace with experimental data!")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_samples', type=int, default=50000)
    parser.add_argument('--n_trajectories', type=int, default=1000)
    parser.add_argument('--trajectory_length', type=int, default=50)
    parser.add_argument('--output_dir', type=str, default='data/')
    parser.add_argument('--config_dim', type=int, default=6)
    parser.add_argument('--state_dim', type=int, default=24)
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate kinematics data
    generate_kinematics_data(
        args.n_samples,
        str(output_dir / 'kinematics_data.npz'),
        args.config_dim
    )
    
    # Generate dynamics data
    generate_dynamics_data(
        args.n_trajectories,
        args.trajectory_length,
        str(output_dir / 'dynamics_data.npz'),
        args.state_dim
    )


if __name__ == '__main__':
    main()
