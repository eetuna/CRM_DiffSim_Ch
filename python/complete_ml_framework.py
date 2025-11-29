
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from pathlib import Path
import pickle
from typing import List, Dict, Optional
import sys
sys.path.insert(0, '../build_debug')

import crm_cpp

# ============================================================================
# PART A: KINEMATICS MODELS
# ============================================================================

# ------------------------------
# A1. FULL ML KINEMATICS
# ------------------------------

class FullMLKinematicsDataset(Dataset):
    """Dataset for pure ML kinematics: currents → tip position (no physics model)."""
    
    def __init__(self, currents: np.ndarray, tip_positions: np.ndarray,
                 normalize: bool = True):
        self.currents = currents.copy()
        self.tip_positions = tip_positions.copy()
        
        # Normalization stats
        self.currents_mean = currents.mean(axis=0)
        self.currents_std = currents.std(axis=0) + 1e-8
        self.tip_mean = tip_positions.mean(axis=0)
        self.tip_std = tip_positions.std(axis=0) + 1e-8
        
        if normalize:
            self.currents = (self.currents - self.currents_mean) / self.currents_std
            self.tip_positions = (self.tip_positions - self.tip_mean) / self.tip_std
        
        self.currents = torch.from_numpy(self.currents).float()
        self.tip_positions = torch.from_numpy(self.tip_positions).float()
    
    def __len__(self):
        return len(self.currents)
    
    def __getitem__(self, idx):
        return {
            'currents': self.currents[idx],
            'tip_position': self.tip_positions[idx]
        }


class FullMLKinematics(nn.Module):
    """
    Pure ML forward kinematics: currents → tip position.
    No physics model involved.
    """
    
    def __init__(self, input_dim: int = 3, hidden_dim: int = 128, output_dim: int = 3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
    
    def forward(self, currents):
        return self.net(currents)


# ------------------------------
# A2. RESIDUAL KINEMATICS (Physics + ML)
# ------------------------------

class ResidualKinematicsDataset(Dataset):
    """Dataset for residual kinematics: (currents, model_pred) → residual."""
    
    def __init__(self, currents: np.ndarray, model_predictions: np.ndarray,
                 true_positions: np.ndarray, normalize: bool = True):
        self.currents = currents.copy()
        self.model_pred = model_predictions.copy()
        self.residuals = true_positions - model_predictions  # What physics missed
        
        # Normalization stats
        self.currents_mean = currents.mean(axis=0)
        self.currents_std = currents.std(axis=0) + 1e-8
        self.residual_mean = self.residuals.mean(axis=0)
        self.residual_std = self.residuals.std(axis=0) + 1e-8
        
        if normalize:
            self.currents = (self.currents - self.currents_mean) / self.currents_std
            self.residuals = (self.residuals - self.residual_mean) / self.residual_std
        
        self.currents = torch.from_numpy(self.currents).float()
        self.model_pred = torch.from_numpy(self.model_pred).float()
        self.residuals = torch.from_numpy(self.residuals).float()
    
    def __len__(self):
        return len(self.currents)
    
    def __getitem__(self, idx):
        return {
            'currents': self.currents[idx],
            'model_pred': self.model_pred[idx],
            'residual': self.residuals[idx]
        }


class ResidualKinematics(nn.Module):
    """
    Residual kinematics: learns correction to physics model.
    Final prediction = physics_pred + NN(currents, physics_pred)
    """
    
    def __init__(self, input_dim: int = 6, hidden_dim: int = 64, output_dim: int = 3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
    
    def forward(self, currents, model_pred):
        x = torch.cat([currents, model_pred], dim=-1)
        return self.net(x)


# ============================================================================
# PART B: DYNAMICS MODELS
# ============================================================================

# ------------------------------
# B1. FULL ML DYNAMICS
# ------------------------------

class FullMLDynamicsDataset(Dataset):
    """
    Dataset for pure ML dynamics: (state_t, action_t) → state_{t+1}.
    No physics model involved.
    """
    
    def __init__(self, trajectories: List[Dict], seq_length: int = 5,
                 normalize: bool = True):
        self.data = []
        
        # Collect all states and actions for normalization
        all_states = np.concatenate([t['states'] for t in trajectories], axis=0)
        all_currents = np.concatenate([t['currents'] for t in trajectories], axis=0)
        
        self.state_mean = all_states.mean(axis=0)
        self.state_std = all_states.std(axis=0) + 1e-8
        self.currents_mean = all_currents.mean(axis=0)
        self.currents_std = all_currents.std(axis=0) + 1e-8
        
        # Build dataset
        for traj in trajectories:
            states = traj['states']
            currents = traj['currents']
            
            if normalize:
                states = (states - self.state_mean) / self.state_std
                currents = (currents - self.currents_mean) / self.currents_std
            
            for t in range(seq_length, len(states) - 1):
                self.data.append({
                    'state_history': states[t-seq_length:t].copy(),
                    'current_history': currents[t-seq_length:t].copy(),
                    'state_t': states[t].copy(),
                    'current_t': currents[t].copy(),
                    'state_next': states[t+1].copy()
                })
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        item = self.data[idx]
        return {k: torch.from_numpy(v).float() for k, v in item.items()}


class FullMLDynamicsLSTM(nn.Module):
    """
    Pure ML dynamics using LSTM: (state_history, action_history) → next_state.
    No physics model involved.
    """
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3,
                 hidden_dim: int = 64, num_layers: int = 2):
        super().__init__()
        self.lstm = nn.LSTM(
            input_size=state_dim + action_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True
        )
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_history, current_history, state_t, current_t):
        combined = torch.cat([state_history, current_history], dim=-1)
        _, (h_n, _) = self.lstm(combined)
        return self.fc(h_n[-1])


class FullMLDynamicsMLP(nn.Module):
    """
    Pure ML dynamics using MLP: (state_t, action_t) → state_{t+1}.
    Ignores history (memoryless assumption).
    """
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3, hidden_dim: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim + action_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_history, current_history, state_t, current_t):
        x = torch.cat([state_t, current_t], dim=-1)
        return self.net(x)


# ------------------------------
# B2. RESIDUAL DYNAMICS (Physics + ML)
# ------------------------------

class ResidualDynamicsDataset(Dataset):
    """
    Dataset for residual dynamics: 
    (state_t, action_t, physics_pred_{t+1}) → residual.
    """
    
    def __init__(self, trajectories: List[Dict], 
                 physics_predictions: List[np.ndarray],
                 seq_length: int = 5, normalize: bool = True):
        """
        Args:
            trajectories: list of trajectory dicts with 'states', 'currents'
            physics_predictions: list of (T, 3) arrays with physics model predictions
        """
        self.data = []
        
        # Normalization stats
        all_states = np.concatenate([t['states'] for t in trajectories], axis=0)
        all_currents = np.concatenate([t['currents'] for t in trajectories], axis=0)
        all_residuals = []
        
        for traj, phys_pred in zip(trajectories, physics_predictions):
            residual = traj['states'][1:] - phys_pred[:-1]  # true - physics
            all_residuals.append(residual)
        
        all_residuals = np.concatenate(all_residuals, axis=0)
        
        self.state_mean = all_states.mean(axis=0)
        self.state_std = all_states.std(axis=0) + 1e-8
        self.currents_mean = all_currents.mean(axis=0)
        self.currents_std = all_currents.std(axis=0) + 1e-8
        self.residual_mean = all_residuals.mean(axis=0)
        self.residual_std = all_residuals.std(axis=0) + 1e-8
        
        # Build dataset
        for traj, phys_pred in zip(trajectories, physics_predictions):
            states = traj['states']
            currents = traj['currents']
            
            if normalize:
                states_norm = (states - self.state_mean) / self.state_std
                currents_norm = (currents - self.currents_mean) / self.currents_std
            else:
                states_norm = states
                currents_norm = currents
            
            for t in range(seq_length, len(states) - 1):
                physics_next = phys_pred[t]  # Physics prediction for t+1
                true_next = states[t+1]
                residual = true_next - physics_next
                
                if normalize:
                    residual = (residual - self.residual_mean) / self.residual_std
                
                self.data.append({
                    'state_history': states_norm[t-seq_length:t].copy(),
                    'current_history': currents_norm[t-seq_length:t].copy(),
                    'state_t': states_norm[t].copy(),
                    'current_t': currents_norm[t].copy(),
                    'physics_pred': physics_next.copy(),
                    'residual': residual.copy()
                })
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        item = self.data[idx]
        return {k: torch.from_numpy(v).float() for k, v in item.items()}


class ResidualDynamicsMLP(nn.Module):
    """
    Residual dynamics (MLP): learns correction to physics model.
    Final prediction = physics_pred + NN(state_t, action_t, physics_pred)
    
    For your magnetic catheter (no hysteresis), this is likely the best choice.
    """
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3, hidden_dim: int = 64):
        super().__init__()
        # Input: state_t (3) + action_t (3) + physics_pred (3) = 9
        self.net = nn.Sequential(
            nn.Linear(state_dim + action_dim + state_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_t, current_t, physics_pred):
        x = torch.cat([state_t, current_t, physics_pred], dim=-1)
        return self.net(x)


class ResidualDynamicsLSTM(nn.Module):
    """
    Residual dynamics (LSTM): learns correction with history.
    Use this only if you observe time-dependent effects.
    """
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3,
                 hidden_dim: int = 64, num_layers: int = 2):
        super().__init__()
        self.lstm = nn.LSTM(
            input_size=state_dim + action_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True
        )
        # Additional input: physics_pred
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim + state_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_history, current_history, state_t, current_t, physics_pred):
        combined = torch.cat([state_history, current_history], dim=-1)
        _, (h_n, _) = self.lstm(combined)
        h_final = h_n[-1]
        x = torch.cat([h_final, physics_pred], dim=-1)
        return self.fc(x)


# ============================================================================
# PART C: UNIFIED TRAINER
# ============================================================================

class UnifiedTrainer:
    """Trainer that handles all model types."""
    
    def __init__(self, model: nn.Module, model_type: str,
                 device: str = 'cpu', learning_rate: float = 1e-3):
        """
        Args:
            model: PyTorch model
            model_type: one of ['full_kin', 'residual_kin', 'full_dyn', 'residual_dyn']
        """
        self.model = model.to(device)
        self.model_type = model_type
        self.device = device
        self.optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        self.criterion = nn.MSELoss()
        self.train_losses = []
        self.val_losses = []
    
    def _forward(self, batch):
        """Forward pass based on model type."""
        if self.model_type == 'full_kin':
            pred = self.model(batch['currents'])
            target = batch['tip_position']
        
        elif self.model_type == 'residual_kin':
            pred = self.model(batch['currents'], batch['model_pred'])
            target = batch['residual']
        
        elif self.model_type == 'full_dyn':
            pred = self.model(
                batch['state_history'], batch['current_history'],
                batch['state_t'], batch['current_t']
            )
            target = batch['state_next']
        
        elif self.model_type == 'residual_dyn':
            pred = self.model(
                batch['state_t'], batch['current_t'], batch['physics_pred']
            )
            target = batch['residual']
        
        else:
            raise ValueError(f"Unknown model_type: {self.model_type}")
        
        return pred, target
    
    def train_epoch(self, loader: DataLoader) -> float:
        self.model.train()
        total_loss = 0.0
        
        for batch in loader:
            batch = {k: v.to(self.device) for k, v in batch.items()}
            pred, target = self._forward(batch)
            loss = self.criterion(pred, target)
            
            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
            self.optimizer.step()
            
            total_loss += loss.item()
        
        avg_loss = total_loss / len(loader)
        self.train_losses.append(avg_loss)
        return avg_loss
    
    def validate(self, loader: DataLoader) -> float:
        self.model.eval()
        total_loss = 0.0
        
        with torch.no_grad():
            for batch in loader:
                batch = {k: v.to(self.device) for k, v in batch.items()}
                pred, target = self._forward(batch)
                loss = self.criterion(pred, target)
                total_loss += loss.item()
        
        avg_loss = total_loss / len(loader)
        self.val_losses.append(avg_loss)
        return avg_loss
    
    def train(self, train_loader, val_loader, epochs: int = 50, patience: int = 10):
        best_val = float('inf')
        patience_count = 0
        
        for epoch in range(epochs):
            train_loss = self.train_epoch(train_loader)
            val_loss = self.validate(val_loader)
            
            print(f"Epoch {epoch+1}/{epochs} | Train: {train_loss:.6f} | Val: {val_loss:.6f}")
            
            if val_loss < best_val:
                best_val = val_loss
                patience_count = 0
                torch.save(self.model.state_dict(), f"best_{self.model_type}.pth")
            else:
                patience_count += 1
                if patience_count >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break
    
    def plot_loss(self, save_path: str = None):
        plt.figure(figsize=(10, 5))
        plt.plot(self.train_losses, label='Train')
        plt.plot(self.val_losses, label='Val')
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        plt.title(f'{self.model_type} Training')
        plt.legend()
        plt.grid()
        if save_path:
            plt.savefig(save_path, dpi=150)
        plt.show()


# ============================================================================
# PART D: DATA COLLECTION
# ============================================================================

class DataCollector:
    """Unified data collector for all model types."""
    
    def __init__(self, crm_kin: crm_cpp.CosseratRod,
                 crm_dyn: crm_cpp.CosseratRodDynamics,
                 save_dir: str = "./data"):
        self.crm_kin = crm_kin
        self.crm_dyn = crm_dyn
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(exist_ok=True)
    
    def collect_kinematics_data(self, num_samples: int = 500,
                                camera_noise_std: float = 0.5e-3) -> Dict:
        """Collect static kinematics data."""
        print(f"Collecting {num_samples} kinematics samples...")
        
        data = {
            'currents': [],
            'inserted_length': [],
            'model_pred': [],
            'camera_measured': []
        }
        
        for i in range(num_samples):
            currents = np.random.uniform(-0.1, 0.1, 3)
            inserted_length = np.random.uniform(0.07, 0.10)
            
            # Physics model prediction
            model_pred = np.array(
                self.crm_kin.forward_kinematics(currents.tolist(), inserted_length)
            )
            
            # Simulated camera measurement (replace with real data)
            camera_measured = model_pred + np.random.randn(3) * camera_noise_std
            
            data['currents'].append(currents)
            data['inserted_length'].append(inserted_length)
            data['model_pred'].append(model_pred)
            data['camera_measured'].append(camera_measured)
        
        return {k: np.array(v) for k, v in data.items()}
    
    def collect_dynamics_data(self, num_trajectories: int = 50,
                              traj_length: int = 100,
                              dt: float = 0.01) -> List[Dict]:
        """Collect dynamics trajectory data with physics predictions."""
        print(f"Collecting {num_trajectories} trajectories...")
        
        trajectories = []
        physics_predictions = []
        
        for traj_idx in range(num_trajectories):
            self.crm_dyn.reset()
            
            states = []
            currents = []
            phys_preds = []
            
            for step in range(traj_length):
                current = np.random.uniform(-0.05, 0.05, 3)
                
                # Get current state before stepping
                state_before = np.array(self.crm_dyn.get_state()[:3])
                
                # Physics model prediction for next state
                # (In practice, you'd run a single-step physics prediction here)
                phys_pred = state_before + np.random.randn(3) * 0.001  # Placeholder
                
                # Actual dynamics step
                self.crm_dyn.step(current.tolist(), dt)
                state_after = np.array(self.crm_dyn.get_state()[:3])
                
                states.append(state_before)
                currents.append(current)
                phys_preds.append(phys_pred)
            
            trajectories.append({
                'states': np.array(states),
                'currents': np.array(currents),
                'dt': dt
            })
            physics_predictions.append(np.array(phys_preds))
        
        return trajectories, physics_predictions
    
    def save_data(self, data, filename: str):
        with open(self.save_dir / filename, 'wb') as f:
            pickle.dump(data, f)
        print(f"Saved to {self.save_dir / filename}")


# ============================================================================
# PART E: MAIN TRAINING SCRIPTS
# ============================================================================

def train_all_kinematics_models(kin_data: Dict, device: str = 'cpu'):
    """Train both full ML and residual kinematics models."""
    
    print("\n" + "="*60)
    print("TRAINING KINEMATICS MODELS")
    print("="*60)
    
    split_idx = int(0.8 * len(kin_data['currents']))
    
    # ----- Full ML Kinematics -----
    print("\n--- Full ML Kinematics ---")
    full_dataset = FullMLKinematicsDataset(
        kin_data['currents'],
        kin_data['camera_measured']
    )
    
    train_loader = DataLoader(
        torch.utils.data.Subset(full_dataset, range(split_idx)),
        batch_size=32, shuffle=True
    )
    val_loader = DataLoader(
        torch.utils.data.Subset(full_dataset, range(split_idx, len(full_dataset))),
        batch_size=32, shuffle=False
    )
    
    full_model = FullMLKinematics(input_dim=3, hidden_dim=128, output_dim=3)
    trainer = UnifiedTrainer(full_model, 'full_kin', device=device)
    trainer.train(train_loader, val_loader, epochs=50)
    trainer.plot_loss('full_kin_loss.png')
    
    # ----- Residual Kinematics -----
    print("\n--- Residual Kinematics ---")
    res_dataset = ResidualKinematicsDataset(
        kin_data['currents'],
        kin_data['model_pred'],
        kin_data['camera_measured']
    )
    
    train_loader = DataLoader(
        torch.utils.data.Subset(res_dataset, range(split_idx)),
        batch_size=32, shuffle=True
    )
    val_loader = DataLoader(
        torch.utils.data.Subset(res_dataset, range(split_idx, len(res_dataset))),
        batch_size=32, shuffle=False
    )
    
    res_model = ResidualKinematics(input_dim=6, hidden_dim=64, output_dim=3)
    trainer = UnifiedTrainer(res_model, 'residual_kin', device=device)
    trainer.train(train_loader, val_loader, epochs=50)
    trainer.plot_loss('residual_kin_loss.png')
    
    return {
        'full_ml': full_model,
        'residual': res_model,
        'full_dataset': full_dataset,
        'res_dataset': res_dataset
    }


def train_all_dynamics_models(trajectories: List[Dict], 
                              physics_predictions: List[np.ndarray],
                              device: str = 'cpu'):
    """Train all dynamics models: Full ML (MLP, LSTM) and Residual (MLP, LSTM)."""
    
    print("\n" + "="*60)
    print("TRAINING DYNAMICS MODELS")
    print("="*60)
    
    # ----- Full ML Dynamics (MLP) -----
    print("\n--- Full ML Dynamics (MLP) ---")
    full_dataset = FullMLDynamicsDataset(trajectories, seq_length=5)
    
    split_idx = int(0.8 * len(full_dataset))
    train_loader = DataLoader(
        torch.utils.data.Subset(full_dataset, range(split_idx)),
        batch_size=32, shuffle=True
    )
    val_loader = DataLoader(
        torch.utils.data.Subset(full_dataset, range(split_idx, len(full_dataset))),
        batch_size=32, shuffle=False
    )
    
    full_mlp = FullMLDynamicsMLP(state_dim=3, action_dim=3, hidden_dim=128)
    trainer = UnifiedTrainer(full_mlp, 'full_dyn', device=device)
    trainer.train(train_loader, val_loader, epochs=50)
    trainer.plot_loss('full_dyn_mlp_loss.png')
    
    # ----- Full ML Dynamics (LSTM) -----
    print("\n--- Full ML Dynamics (LSTM) ---")
    full_lstm = FullMLDynamicsLSTM(state_dim=3, action_dim=3, hidden_dim=64)
    trainer = UnifiedTrainer(full_lstm, 'full_dyn', device=device)
    trainer.train(train_loader, val_loader, epochs=50)
    trainer.plot_loss('full_dyn_lstm_loss.png')
    
    # ----- Residual Dynamics (MLP) -----
    print("\n--- Residual Dynamics (MLP) ---")
    res_dataset = ResidualDynamicsDataset(trajectories, physics_predictions, seq_length=5)
    
    split_idx = int(0.8 * len(res_dataset))
    train_loader = DataLoader(
        torch.utils.data.Subset(res_dataset, range(split_idx)),
        batch_size=32, shuffle=True
    )
    val_loader = DataLoader(
        torch.utils.data.Subset(res_dataset, range(split_idx, len(res_dataset))),
        batch_size=32, shuffle=False
    )
    
    # Custom trainer for residual dynamics
    res_mlp = ResidualDynamicsMLP(state_dim=3, action_dim=3, hidden_dim=64)
    res_mlp = res_mlp.to(device)
    optimizer = optim.Adam(res_mlp.parameters(), lr=1e-3)
    criterion = nn.MSELoss()
    
    print("Training Residual Dynamics MLP...")
    for epoch in range(50):
        res_mlp.train()
        train_loss = 0.0
        for batch in train_loader:
            batch = {k: v.to(device) for k, v in batch.items()}
            pred = res_mlp(batch['state_t'], batch['current_t'], batch['physics_pred'])
            loss = criterion(pred, batch['residual'])
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}: Loss = {train_loss/len(train_loader):.6f}")
    
    torch.save(res_mlp.state_dict(), 'best_residual_dyn.pth')
    
    return {
        'full_mlp': full_mlp,
        'full_lstm': full_lstm,
        'residual_mlp': res_mlp,
        'datasets': {
            'full': full_dataset,
            'residual': res_dataset
        }
    }


# ============================================================================
# PART F: MAIN ENTRY POINT
# ============================================================================

def main():
    """Train all 4 model types."""
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")
    
    # Initialize C++ models
    crm_kin = crm_cpp.CosseratRod("config.json")
    crm_dyn = crm_cpp.CosseratRodDynamics("config.json")
    
    collector = DataCollector(crm_kin, crm_dyn)
    
    # Collect data
    kin_data = collector.collect_kinematics_data(num_samples=500)
    trajectories, physics_preds = collector.collect_dynamics_data(num_trajectories=50)
    
    # Train kinematics models
    kin_models = train_all_kinematics_models(kin_data, device=device)
    
    # Train dynamics models
    dyn_models = train_all_dynamics_models(trajectories, physics_preds, device=device)
    
    print("\n" + "="*60)
    print("ALL MODELS TRAINED SUCCESSFULLY")
    print("="*60)
    print("\nSaved models:")
    print("  - best_full_kin.pth (Full ML Kinematics)")
    print("  - best_residual_kin.pth (Residual Kinematics)")
    print("  - best_full_dyn.pth (Full ML Dynamics)")
    print("  - best_residual_dyn.pth (Residual Dynamics)")


if __name__ == "__main__":
    main()
