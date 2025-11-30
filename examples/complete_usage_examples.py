"""
Complete usage examples for CRM Catheter Python API
"""
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from crm_wrapper import CatheterModel
from envs.catheter_env import CatheterEnv
from models.catheter_dataset import CatheterDataset, CatheterResNet, train_forward_model
from models.contact_models import ContactForceEstimator, ResidualContactDynamics, ContactDataset, train_contact_force_model
from rl.td3_sac import TD3Agent, SACAgent, ReplayBuffer, train_rl_agent
from rl.sim_to_real import DomainRandomization, ModelBasedSimulator, collect_real_world_data

def example_forward_kinematics():
    """Example: Basic forward kinematics"""
    print("=" * 60)
    print("Forward Kinematics Example")
    print("=" * 60)
    
    catheter = CatheterModel()
    
    currents = np.array([[1.0, 0.5, 0.2]])  # 1 actuator, 3 coils
    insertion_length = 80.0  # mm
    
    result = catheter.forward_kinematics(currents, insertion_length)
    
    print(f"Tip Position: {result['tip_position']}")
    print(f"Tip Rotation:\n{result['tip_rotation']}")
    print(f"Converged: {result['converged']}")
    print(f"Energy: {result['potential_energy']:.4f}")
    print()

def example_jacobian():
    """Example: Jacobian computation for control"""
    print("=" * 60)
    print("Jacobian Example")
    print("=" * 60)
    
    catheter = CatheterModel()
    
    currents = np.array([[1.0, 0.5, 0.2]])
    insertion_length = 80.0
    
    J = catheter.compute_jacobian(currents, insertion_length, analytical=True)
    
    print(f"Jacobian shape: {J.shape}")
    print(f"Jacobian:\n{J}")
    print()

def example_rl_environment():
    """Example: RL environment for catheter control"""
    print("=" * 60)
    print("RL Environment Example")
    print("=" * 60)
    
    env = CatheterEnv(max_steps=200, target_radius=2.0)
    
    state, _ = env.reset()
    print(f"Initial state shape: {state.shape}")
    
    for step in range(10):
        action = env.action_space.sample()
        next_state, reward, terminated, truncated, info = env.step(action)
        
        if step == 0:
            print(f"Action: {action}")
            print(f"Reward: {reward:.4f}")
            print(f"Next state: {next_state[:5]}...")  # First 5 elements
        
        if terminated or truncated:
            break
    
    print()

def example_batch_processing():
    """Example: Efficient batch processing for ML"""
    print("=" * 60)
    print("Batch Processing Example")
    print("=" * 60)
    
    catheter = CatheterModel()
    
    batch_size = 1000
    currents_batch = np.random.uniform(-2, 2, (batch_size, 3))
    insertions_batch = np.random.uniform(10, 110, (batch_size,))
    
    print(f"Processing {batch_size} samples...")
    outputs = catheter.batch_forward_kinematics(currents_batch, insertions_batch)
    
    print(f"Output shape: {outputs.shape}")
    print(f"Sample output: {outputs[0]}")
    print()

def example_supervised_learning():
    """Example: Train forward kinematics neural network"""
    print("=" * 60)
    print("Supervised Learning Example")
    print("=" * 60)
    
    batch_size = 64
    n_samples = 1000
    
    inputs = torch.randn(n_samples, 4)  # currents + insertion
    targets = torch.randn(n_samples, 15)  # position + rotation + curvature
    
    dataset = CatheterDataset(inputs, targets)
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    
    model = CatheterResNet(input_dim=4, output_dim=15)
    
    print("Training forward kinematics model...")
    trained_model = train_forward_model(model, train_loader, val_loader, n_epochs=5)
    print("Model trained successfully")
    print()

def example_contact_force_estimation():
    """Example: Train contact force estimator"""
    print("=" * 60)
    print("Contact Force Estimation Example")
    print("=" * 60)
    
    contact_model = ContactForceEstimator(config_dim=6, actuation_dim=4, env_dim=8)
    
    n_trajectories = 50
    seq_len = 100
    trajectories = []
    
    for _ in range(n_trajectories):
        states = np.random.randn(seq_len, 6)
        controls = np.random.randn(seq_len, 4)
        contact_forces = np.random.randn(seq_len, 3)
        env_features = np.random.randn(seq_len, 8)
        in_contact = (np.random.rand(seq_len) > 0.7).astype(float)
        contact_labels = in_contact.astype(float)
        
        trajectories.append({
            'states': states,
            'controls': controls,
            'contact_forces': contact_forces,
            'env_features': env_features,
            'contact_labels': contact_labels
        })
    
    dataset = ContactDataset(trajectories)
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])
    
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=32)
    
    print("Training contact force model...")
    trained_model = train_contact_force_model(contact_model, train_loader, val_loader, n_epochs=10)
    print("Contact force model trained successfully")
    print()

def example_td3_training():
    """Example: Train TD3 agent for catheter control"""
    print("=" * 60)
    print("TD3 Training Example")
    print("=" * 60)
    
    env = CatheterEnv(max_steps=200, target_radius=2.0)
    
    state_dim = env.observation_space.shape[0]
    action_dim = env.action_space.shape[0]
    
    agent = TD3Agent(
        state_dim=state_dim,
        action_dim=action_dim,
        hidden_dim=256,
        lr=3e-4,
        gamma=0.99,
        tau=0.005,
        policy_noise=0.2,
        noise_clip=0.5,
        policy_freq=2
    )
    
    replay_buffer = ReplayBuffer(capacity=100000)
    
    print("Training TD3 agent...")
    results = train_rl_agent(
        env=env,
        agent=agent,
        replay_buffer=replay_buffer,
        n_episodes=20,  # Use 1000 for real training
        max_steps=200,
        batch_size=256,
        warmup_steps=5000,
        eval_freq=10,
        save_freq=50
    )
    
    print(f"TD3 training completed")
    print(f"Final eval reward: {results['eval_rewards'][-1]:.2f}")
    print()

def example_sac_training():
    """Example: Train SAC agent for catheter control"""
    print("=" * 60)
    print("SAC Training Example")
    print("=" * 60)
    
    env = CatheterEnv(max_steps=200, target_radius=2.0)
    
    state_dim = env.observation_space.shape[0]
    action_dim = env.action_space.shape[0]
    
    agent = SACAgent(
        state_dim=state_dim,
        action_dim=action_dim,
        hidden_dim=256,
        lr=3e-4,
        gamma=0.99,
        tau=0.005,
        alpha=0.2,
        auto_tune_alpha=True
    )
    
    replay_buffer = ReplayBuffer(capacity=100000)
    
    print("Training SAC agent with automatic temperature tuning...")
    results = train_rl_agent(
        env=env,
        agent=agent,
        replay_buffer=replay_buffer,
        n_episodes=20,
        max_steps=200,
        batch_size=256,
        warmup_steps=5000,
        eval_freq=10,
        save_freq=50
    )
    
    print(f"SAC training completed")
    print(f"Final eval reward: {results['eval_rewards'][-1]:.2f}")
    print(f"Final alpha: {agent.alpha:.4f}")
    print()

def example_domain_randomization():
    """Example: Sim-to-real with domain randomization"""
    print("=" * 60)
    print("Domain Randomization Example")
    print("=" * 60)
    
    param_ranges = {
        'mass': (0.01, 0.05),
        'damping': (0.01, 0.1),
        'stiffness': (0.8, 1.2),
        'friction': (0.1, 0.5)
    }
    
    randomizer = DomainRandomization(param_ranges)
    
    print("Generating randomized environments...")
    for i in range(5):
        params = randomizer.sample_parameters()
        print(f"Env {i+1}: {params}")
    
    print()

if __name__ == "__main__":
    example_forward_kinematics()
    example_jacobian()
    example_rl_environment()
    example_batch_processing()
    example_supervised_learning()
    example_contact_force_estimation()
    example_td3_training()
    example_sac_training()
    example_domain_randomization()
