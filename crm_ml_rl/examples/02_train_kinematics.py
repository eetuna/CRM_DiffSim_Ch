"""
Train kinematics residual model.

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

from crm_ml_rl.models import KinematicsResidualNN


def train_kinematics_model(
    data_path: str,
    epochs: int,
    batch_size: int,
    lr: float,
    save_path: str
):
    """Train kinematics residual model."""
    
    # Load data
    print(f"Loading data from {data_path}...")
    data = np.load(data_path)
    
    angles = data['angles']
    currents = data['currents']
    positions = data['positions']
    
    print(f"Loaded {len(angles)} samples")
    print(f"Config dim: {angles.shape[1]}, Currents: {currents.shape[1]}")
    
    # For now, simulate residuals as zero (no real data)
    # In practice, residuals = positions_real - positions_analytical
    residuals = np.zeros_like(positions)
    
    print("\n" + "="*60)
    print("WARNING: Training on ZERO residuals (no real data)")
    print("For actual residual learning, you need experimental measurements!")
    print("="*60 + "\n")
    
    # Prepare data
    X = np.concatenate([currents, angles], axis=1)
    y = residuals
    
    input_dim = X.shape[1]
    print(f"Input dimension: {input_dim}")
    
    # Split train/val
    split = int(0.9 * len(X))
    X_train, X_val = X[:split], X[split:]
    y_train, y_val = y[:split], y[split:]
    
    # Create dataloaders
    train_dataset = TensorDataset(
        torch.FloatTensor(X_train),
        torch.FloatTensor(y_train)
    )
    val_dataset = TensorDataset(
        torch.FloatTensor(X_val),
        torch.FloatTensor(y_val)
    )
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size)
    
    # Model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    model = KinematicsResidualNN(input_dim=input_dim).to(device)
    
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()
    
    # Training loop
    print(f"\nTraining for {epochs} epochs...")
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        
        for batch_x, batch_y in tqdm(train_loader, desc=f'Epoch {epoch+1}/{epochs}'):
            batch_x = batch_x.to(device)
            batch_y = batch_y.to(device)
            
            optimizer.zero_grad()
            pred = model(batch_x)
            loss = criterion(pred, batch_y)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch_x, batch_y in val_loader:
                batch_x = batch_x.to(device)
                batch_y = batch_y.to(device)
                pred = model(batch_x)
                val_loss += criterion(pred, batch_y).item()
        
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
    parser.add_argument('--data_path', type=str, default='data/kinematics_data.npz')
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch_size', type=int, default=512)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--save_path', type=str, default='models/kinematics_residual.pth')
    args = parser.parse_args()
    
    train_kinematics_model(
        args.data_path,
        args.epochs,
        args.batch_size,
        args.lr,
        args.save_path
    )


if __name__ == '__main__':
    main()
