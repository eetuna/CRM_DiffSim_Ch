"""Data generation utilities."""

import numpy as np
import h5py
from pathlib import Path
from tqdm import tqdm

from .kinematics import ForwardKinematicsAPI
from .dynamics import DynamicsAPI


class KinematicsDataGenerator:
    """Generate forward kinematics training data."""
    
    def __init__(self, n_segments: int = 3):
        self.fk_api = ForwardKinematicsAPI(n_segments=n_segments)
        self.n_dof = n_segments * 3
        
    def generate(self, n_samples: int, output_file: str, 
                 q_range: tuple = (-np.pi/2, np.pi/2),
                 I_range: tuple = (-2.0, 2.0)):
        """
        Generate FK training data.
        
        Args:
            n_samples: Number of samples to generate
            output_file: Output HDF5 file path
            q_range: Joint angle range (min, max)
            I_range: Current range (min, max)
        """
        print(f"\nGenerating {n_samples} kinematics samples...")
        
        q_data = np.zeros((n_samples, self.n_dof))
        I_data = np.zeros((n_samples, 6))
        tip_positions = np.zeros((n_samples, 3))
        
        for i in tqdm(range(n_samples)):
            q = np.random.uniform(q_range[0], q_range[1], size=self.n_dof)
            I = np.random.uniform(I_range[0], I_range[1], size=6)
            
            try:
                tip_pos = self.fk_api.get_tip_position(q, I)
                q_data[i] = q
                I_data[i] = I
                tip_positions[i] = tip_pos
            except Exception as e:
                print(f"Warning: Sample {i} failed: {e}")
                continue
        
        # Save to HDF5
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_file, 'w') as f:
            f.create_dataset('q', data=q_data, compression='gzip')
            f.create_dataset('I', data=I_data, compression='gzip')
            f.create_dataset('tip_positions', data=tip_positions, compression='gzip')
            f.attrs['n_samples'] = n_samples
            f.attrs['n_dof'] = self.n_dof
            f.attrs['q_range'] = q_range
            f.attrs['I_range'] = I_range
        
        print(f"✓ Saved {n_samples} samples to {output_file}")
        return output_file


class DynamicsDataGenerator:
    """Generate dynamics training data."""
    
    def __init__(self, state_dim: int = 18, control_dim: int = 6, dt: float = 0.001):
        self.dyn_api = DynamicsAPI(state_dim=state_dim, control_dim=control_dim, dt=dt)
        self.state_dim = state_dim
        self.control_dim = control_dim
        
    def generate(self, n_trajectories: int, steps_per_traj: int, 
                 output_file: str,
                 control_range: tuple = (-2.0, 2.0),
                 initial_state_std: float = 0.1):
        """
        Generate dynamics training data.
        
        Args:
            n_trajectories: Number of trajectories
            steps_per_traj: Steps per trajectory
            output_file: Output HDF5 file path
            control_range: Control input range
            initial_state_std: Initial state standard deviation
        """
        print(f"\nGenerating {n_trajectories} dynamics trajectories...")
        
        states_list = []
        controls_list = []
        next_states_list = []
        
        for _ in tqdm(range(n_trajectories)):
            # Random initial state
            state = np.random.randn(self.state_dim) * initial_state_std
            
            for _ in range(steps_per_traj):
                # Random control
                control = np.random.uniform(control_range[0], control_range[1], size=self.control_dim)
                
                try:
                    next_state = self.dyn_api.step(state, control)
                    
                    states_list.append(state)
                    controls_list.append(control)
                    next_states_list.append(next_state)
                    
                    state = next_state
                except Exception as e:
                    print(f"Warning: Step failed: {e}")
                    break
        
        states_data = np.array(states_list)
        controls_data = np.array(controls_list)
        next_states_data = np.array(next_states_list)
        
        # Save to HDF5
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_file, 'w') as f:
            f.create_dataset('states', data=states_data, compression='gzip')
            f.create_dataset('controls', data=controls_data, compression='gzip')
            f.create_dataset('next_states', data=next_states_data, compression='gzip')
            f.attrs['n_trajectories'] = n_trajectories
            f.attrs['steps_per_traj'] = steps_per_traj
            f.attrs['state_dim'] = self.state_dim
            f.attrs['control_dim'] = self.control_dim
            f.attrs['control_range'] = control_range
        
        print(f"✓ Saved {len(states_data)} transitions to {output_file}")
        return output_file


def generate_all_data(kinematics_samples: int = 100000,
                     dynamics_trajectories: int = 1000,
                     dynamics_steps: int = 100,
                     output_dir: str = 'data'):
    """
    Generate all training and validation data.
    
    Args:
        kinematics_samples: Number of FK samples for training
        dynamics_trajectories: Number of dynamics trajectories
        dynamics_steps: Steps per trajectory
        output_dir: Output directory
    """
    print("=" * 60)
    print("Data Generation")
    print("=" * 60)
    
    # Kinematics data
    kin_gen = KinematicsDataGenerator(n_segments=3)
    kin_gen.generate(
        n_samples=kinematics_samples,
        output_file=f'{output_dir}/kinematics_train.h5'
    )
    kin_gen.generate(
        n_samples=kinematics_samples // 10,
        output_file=f'{output_dir}/kinematics_val.h5'
    )
    
    # Dynamics data
    dyn_gen = DynamicsDataGenerator(state_dim=18, control_dim=6, dt=0.001)
    dyn_gen.generate(
        n_trajectories=dynamics_trajectories,
        steps_per_traj=dynamics_steps,
        output_file=f'{output_dir}/dynamics_train.h5'
    )
    dyn_gen.generate(
        n_trajectories=dynamics_trajectories // 10,
        steps_per_traj=dynamics_steps,
        output_file=f'{output_dir}/dynamics_val.h5'
    )
    
    print("\n✓ All data generation complete!")
