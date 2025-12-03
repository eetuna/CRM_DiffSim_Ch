"""
Train dynamics model.

NOTE: This trains on analytical data only.
For real residual learning, you need experimental measurements!
"""

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import argparse
from pathlib import Path
from tqdm import tqdm

from crm_ml_rl.models import DynamicsModelNN


def train_dynamics_model(
    data_path: str,
    epochs: int,
    batch_size: int,
    lr: float,
    save_path: str
):
    """Train dynamics model."""
    
    # Load data
    print(f"Loading data from {data_path}...")
    data = np.load(data_path)
    
    states = data['states']
    actions = data['actions']
    next_states = data['next_states']
    
    print(f"Loaded {len(states)} samples")
    print(f"State dim: {states.shape[1]}, Action dim: {actions.shape[1]}")
    
    # For now, simulate corrections as zero (no real data)
    corrections = np.zeros_like(states)
    
    print("\n" + "="*60)
    print("WARNING: Training on ZERO corrections (no real data)")
    print("For actual residual learning, you need experimental measurements!")
    print("="*60 + "\n")
    
    # Split train/val
    split = int(0.9 * len(states))
    states_train, states_val = states[:split], states[split:]
    actions_train, actions_val = actions[:split], actions[split:]
    corrections_train, corrections_val = corrections[:split], corrections[split:]
    
    # Create dataloaders
    train_dataset = TensorDataset(
        torch.FloatTensor(states_train),
        torch.FloatTensor(actions_train),
        torch.FloatTensor(corrections_train)
    )
    val_dataset = TensorDataset(
        torch.FloatTensor(states_val),
        torch.FloatTensor(actions_val),
        torch.FloatTensor(corrections_val)
    )
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    
    # Model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    state_dim = states.shape[1]
    action_dim = actions.shape[1]
    
    model = DynamicsModelNN(
        state_dim=state_dim,
        action_dim=action_dim
    ).to(device)
    
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()
    
    # Training loop
    print(f"\nTraining for {epochs} epochs...")
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        
        for state, action, correction in tqdm(train_loader, desc=f'Epoch {epoch+1}/{epochs}'):
            state = state.to(device)
            action = action.to(device)
            correction = correction.to(device)
            
            optimizer.zero_grad()
            pred = model(state, action)
            loss = criterion(pred, correction)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for state, action, correction in val_loader:
                state = state.to(device)
                action = action.to(device)
                correction = correction.to(device)
                pred = model(state, action)
                val_loss += criterion(pred, correction).item()
        
        train_loss /= len(train_loader)
        val_loss /= len(val_loader)
        
        print(f"Epoch {epoch+1}: Train Loss = {train_loss:.6f}, Val Loss = {val_loss:.6f}")
        
        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            torch.save(model.state_dict(), save_path)
            print(f"  → Best model saved!")
    
    print(f"\nTraining complete! Best val loss: {best_val_loss:.6f}")
    print(f"Model saved to {save_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_path', type=str, default='data/dynamics_data.npz')
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch_size', type=int, default=256)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--save_path', type=str, default='models/dynamics_model.pth')
    args = parser.parse_args()
    
    train_dynamics_model(
        args.data_path,
        args.epochs,
        args.batch_size,
        args.lr,
        args.save_path
    )


if __name__ == '__main__':
    main()
