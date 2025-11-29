import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import crm_cpp
from typing import List, Dict, Tuple

# ============================================================================
# 1. DYNAMICS DATA COLLECTION
# ============================================================================

class DynamicsDataCollector:
    """Collect time-series trajectory data for dynamics learning."""
    
    def __init__(self, crm_dyn: crm_cpp.CosseratRodDynamics, save_dir: str = "./data"):
        self.crm_dyn = crm_dyn
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(exist_ok=True)
    
    def collect_trajectory(self, currents_schedule: np.ndarray, 
                          dt: float = 0.01, sim_time: float = 1.0) -> Dict:
        """
        Collect one trajectory by applying time-varying currents.
        
        Args:
            currents_schedule: (num_steps, 3) array of currents over time
            dt: time step
            sim_time: total simulation time
        
        Returns:
            trajectory: dict with 'states' (positions), 'currents', 'timestamps'
        """
        self.crm_dyn.reset()
        
        states = []  # tip position (3,)
        velocities = []  # tip velocity (3,)
        currents_applied = []
        timestamps = []
        
        num_steps = len(currents_schedule)
        
        for step in range(num_steps):
            current = currents_schedule[step]
            
            # Step the dynamics
            self.crm_dyn.step(current.tolist(), dt=dt)
            
            # Get state: [position (3), orientation_quat (4), lin_vel (3), ang_vel (3)]
            state_full = self.crm_dyn.get_state()  # Returns list of ~13 elements
            
            tip_pos = np.array(state_full[:3])
            tip_vel = np.array(state_full[7:10])  # Adjust indices based on your API
            
            states.append(tip_pos)
            velocities.append(tip_vel)
            currents_applied.append(current)
            timestamps.append(step * dt)
        
        trajectory = {
            'states': np.array(states),           # (num_steps, 3)
            'velocities': np.array(velocities),   # (num_steps, 3)
            'currents': np.array(currents_applied),  # (num_steps, 3)
            'timestamps': np.array(timestamps),    # (num_steps,)
            'dt': dt
        }
        
        return trajectory
    
    def collect_multiple_trajectories(self, num_trajectories: int = 50,
                                      trajectory_length: int = 100,
                                      dt: float = 0.01,
                                      current_scale: float = 0.05) -> List[Dict]:
        """Collect a dataset of random trajectories."""
        
        print(f"Collecting {num_trajectories} trajectories...")
        trajectories = []
        
        for traj_idx in range(num_trajectories):
            # Random time-varying currents (smooth via low-pass filtering)
            raw_currents = np.random.randn(trajectory_length, 3) * current_scale
            
            # Smooth with moving average for realistic time-varying inputs
            window = 5
            smoothed = np.convolve(raw_currents[:, 0], np.ones(window)/window, mode='same')
            currents_schedule = np.zeros_like(raw_currents)
            
            for ch in range(3):
                currents_schedule[:, ch] = np.convolve(
                    raw_currents[:, ch], np.ones(window)/window, mode='same'
                )
            
            # Clip to valid range
            currents_schedule = np.clip(currents_schedule, -0.1, 0.1)
            
            trajectory = self.collect_trajectory(currents_schedule, dt=dt)
            trajectories.append(trajectory)
            
            if (traj_idx + 1) % 10 == 0:
                print(f"  {traj_idx + 1}/{num_trajectories} trajectories collected")
        
        return trajectories
    
    def save_trajectories(self, trajectories: List[Dict], 
                         filename: str = "dynamics_trajectories.pkl"):
        """Save trajectory data."""
        filepath = self.save_dir / filename
        with open(filepath, 'wb') as f:
            pickle.dump(trajectories, f)
        print(f"Trajectories saved to {filepath}")
    
    def load_trajectories(self, filename: str = "dynamics_trajectories.pkl") -> List[Dict]:
        """Load previously collected trajectories."""
        filepath = self.save_dir / filename
        with open(filepath, 'rb') as f:
            trajectories = pickle.load(f)
        print(f"Trajectories loaded from {filepath}")
        return trajectories


# ============================================================================
# 2. DYNAMICS DATASETS
# ============================================================================

class SequencePredictionDataset(Dataset):
    """
    Dataset for learning state transition: 
    (state_t, currents_t) → state_t+1
    
    Also includes sequence history for LSTM.
    """
    
    def __init__(self, trajectories: List[Dict], seq_history_length: int = 5,
                 normalize: bool = True):
        """
        Args:
            trajectories: list of trajectory dicts
            seq_history_length: how many past states/currents to include
            normalize: whether to normalize states/currents
        """
        self.seq_length = seq_history_length
        self.data = []
        
        # Compute normalization stats
        all_states = np.concatenate([t['states'] for t in trajectories], axis=0)
        all_currents = np.concatenate([t['currents'] for t in trajectories], axis=0)
        
        self.state_mean = all_states.mean(axis=0)
        self.state_std = all_states.std(axis=0) + 1e-8
        self.currents_mean = all_currents.mean(axis=0)
        self.currents_std = all_currents.std(axis=0) + 1e-8
        
        # Build dataset
        for traj in trajectories:
            states = traj['states']  # (T, 3)
            currents = traj['currents']  # (T, 3)
            
            # Normalize
            if normalize:
                states_norm = (states - self.state_mean) / self.state_std
                currents_norm = (currents - self.currents_mean) / self.currents_std
            else:
                states_norm = states
                currents_norm = currents
            
            # Create sequences
            for t in range(seq_history_length, len(states) - 1):
                # History: states and currents from t-seq_length to t
                state_history = states_norm[t-seq_history_length:t]  # (seq_length, 3)
                current_history = currents_norm[t-seq_history_length:t]  # (seq_length, 3)
                
                # Current full state + action
                state_t = states_norm[t]  # (3,)
                current_t = currents_norm[t]  # (3,)
                
                # Next state
                state_next = states_norm[t+1]  # (3,)
                
                self.data.append({
                    'state_history': state_history.copy(),
                    'current_history': current_history.copy(),
                    'state_t': state_t.copy(),
                    'current_t': current_t.copy(),
                    'state_next': state_next.copy(),
                })
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        item = self.data[idx]
        return {
            'state_history': torch.from_numpy(item['state_history']).float(),
            'current_history': torch.from_numpy(item['current_history']).float(),
            'state_t': torch.from_numpy(item['state_t']).float(),
            'current_t': torch.from_numpy(item['current_t']).float(),
            'state_next': torch.from_numpy(item['state_next']).float(),
        }


# ============================================================================
# 3. DYNAMICS MODELS
# ============================================================================

class SimpleRNNDynamics(nn.Module):
    """Simple RNN (Elman net) for dynamics."""
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3, 
                 hidden_dim: int = 64):
        super().__init__()
        self.hidden_dim = hidden_dim
        
        # Embed current state + action
        self.input_embed = nn.Linear(state_dim + action_dim, hidden_dim)
        
        # RNN cell
        self.rnn_cell = nn.RNNCell(hidden_dim, hidden_dim)
        
        # Output: predict next state
        self.output = nn.Linear(hidden_dim, state_dim)
    
    def forward(self, state_history, current_history, state_t, current_t):
        """
        Args:
            state_history: (batch, seq_len, 3)
            current_history: (batch, seq_len, 3)
            state_t: (batch, 3)
            current_t: (batch, 3)
        """
        batch_size = state_history.shape[0]
        h = torch.zeros(batch_size, self.hidden_dim, device=state_history.device)
        
        # Process history
        for t in range(state_history.shape[1]):
            s_h = state_history[:, t]  # (batch, 3)
            a_h = current_history[:, t]  # (batch, 3)
            x = torch.cat([s_h, a_h], dim=-1)
            x = self.input_embed(x)
            h = self.rnn_cell(x, h)
        
        # Process current state + action
        x_curr = torch.cat([state_t, current_t], dim=-1)
        x_curr = self.input_embed(x_curr)
        h = self.rnn_cell(x_curr, h)
        
        # Predict next state
        state_next_pred = self.output(h)
        return state_next_pred


class LSTMDynamics(nn.Module):
    """LSTM for dynamics prediction."""
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3, 
                 hidden_dim: int = 64, num_layers: int = 2):
        super().__init__()
        self.state_dim = state_dim
        self.hidden_dim = hidden_dim
        
        # Combine state and action as input
        self.lstm = nn.LSTM(
            input_size=state_dim + action_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True
        )
        
        # Output layer: predict next state
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_history, current_history, state_t, current_t):
        """
        Args:
            state_history: (batch, seq_len, 3)
            current_history: (batch, seq_len, 3)
            state_t: (batch, 3)
            current_t: (batch, 3)
        """
        batch_size = state_history.shape[0]
        
        # Concatenate state and action histories
        combined_history = torch.cat([state_history, current_history], dim=-1)
        
        # LSTM forward
        _, (h_n, _) = self.lstm(combined_history)
        
        # Use last hidden state
        h_final = h_n[-1]  # (batch, hidden_dim)
        
        # Predict next state
        state_next_pred = self.fc(h_final)
        
        return state_next_pred


class NeuralODEDynamics(nn.Module):
    """Neural ODE residual dynamics: state_next = state_t + dt * f_NN(state_t, action_t)"""
    
    def __init__(self, state_dim: int = 3, action_dim: int = 3, 
                 hidden_dim: int = 64, dt: float = 0.01):
        super().__init__()
        self.dt = dt
        
        self.dynamics_net = nn.Sequential(
            nn.Linear(state_dim + action_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_history, current_history, state_t, current_t):
        """
        Simple version: just use current state + action.
        Could also process history if needed.
        """
        x = torch.cat([state_t, current_t], dim=-1)
        delta = self.dynamics_net(x)
        state_next_pred = state_t + self.dt * delta
        return state_next_pred


# ============================================================================
# 4. TRAINING LOOP FOR DYNAMICS
# ============================================================================

class DynamicsTrainer:
    """Trainer for dynamics models."""
    
    def __init__(self, model: nn.Module, device: str = 'cpu',
                 learning_rate: float = 1e-3, weight_decay: float = 1e-5):
        self.model = model.to(device)
        self.device = device
        self.optimizer = optim.Adam(model.parameters(), lr=learning_rate,
                                     weight_decay=weight_decay)
        self.criterion = nn.MSELoss()
        self.train_loss_history = []
        self.val_loss_history = []
    
    def train_epoch(self, train_loader: DataLoader) -> float:
        """Train one epoch."""
        self.model.train()
        total_loss = 0.0
        
        for batch in train_loader:
            # Move to device
            batch = {k: v.to(self.device) for k, v in batch.items()}
            
            # Forward pass
            pred = self.model(
                batch['state_history'],
                batch['current_history'],
                batch['state_t'],
                batch['current_t']
            )
            
            target = batch['state_next']
            loss = self.criterion(pred, target)
            
            # Backward
            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()
            
            total_loss += loss.item()
        
        avg_loss = total_loss / len(train_loader)
        self.train_loss_history.append(avg_loss)
        return avg_loss
    
    def validate(self, val_loader: DataLoader) -> float:
        """Validate on validation set."""
        self.model.eval()
        total_loss = 0.0
        
        with torch.no_grad():
            for batch in val_loader:
                batch = {k: v.to(self.device) for k, v in batch.items()}
                
                pred = self.model(
                    batch['state_history'],
                    batch['current_history'],
                    batch['state_t'],
                    batch['current_t']
                )
                
                target = batch['state_next']
                loss = self.criterion(pred, target)
                total_loss += loss.item()
        
        avg_loss = total_loss / len(val_loader)
        self.val_loss_history.append(avg_loss)
        return avg_loss
    
    def train(self, train_loader: DataLoader, val_loader: DataLoader,
              num_epochs: int = 50, patience: int = 10):
        """Full training with early stopping."""
        best_val_loss = float('inf')
        patience_count = 0
        
        for epoch in range(num_epochs):
            train_loss = self.train_epoch(train_loader)
            val_loss = self.validate(val_loader)
            
            print(f"Epoch {epoch+1}/{num_epochs} | "
                  f"Train Loss: {train_loss:.6f} | Val Loss: {val_loss:.6f}")
            
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_count = 0
                self.save_checkpoint("best_dynamics_model.pth")
            else:
                patience_count += 1
                if patience_count >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break
    
    def save_checkpoint(self, filepath: str):
        torch.save(self.model.state_dict(), filepath)
        print(f"Model saved to {filepath}")
    
    def load_checkpoint(self, filepath: str):
        self.model.load_state_dict(torch.load(filepath, map_location=self.device))
        print(f"Model loaded from {filepath}")
    
    def plot_loss(self):
        plt.figure(figsize=(10, 5))
        plt.plot(self.train_loss_history, label='Train Loss')
        plt.plot(self.val_loss_history, label='Val Loss')
        plt.xlabel('Epoch')
        plt.ylabel('Loss (MSE)')
        plt.legend()
        plt.grid()
        plt.savefig('dynamics_loss_curve.png', dpi=150)
        print("Loss curve saved to dynamics_loss_curve.png")


# ============================================================================
# 5. EVALUATION & ROLLOUT PREDICTION
# ============================================================================

class DynamicsEvaluator:
    """Evaluate dynamics model on rollout predictions."""
    
    def __init__(self, model: nn.Module, dataset: SequencePredictionDataset,
                 device: str = 'cpu'):
        self.model = model.to(device).eval()
        self.dataset = dataset
        self.device = device
    
    def predict_rollout(self, initial_state: torch.Tensor,
                       current_schedule: torch.Tensor,
                       steps: int) -> np.ndarray:
        """
        Predict trajectory by unrolling predictions.
        
        Args:
            initial_state: (3,) tensor
            current_schedule: (steps, 3) tensor
            steps: number of steps to predict
        
        Returns:
            trajectory: (steps, 3) array of predicted states
        """
        trajectory = [initial_state.cpu().numpy()]
        state = initial_state.unsqueeze(0).to(self.device)
        
        # Initialize history with zeros
        seq_len = self.dataset.seq_length
        state_history = torch.zeros(1, seq_len, 3, device=self.device)
        current_history = torch.zeros(1, seq_len, 3, device=self.device)
        
        with torch.no_grad():
            for step in range(steps):
                current = current_schedule[step].unsqueeze(0).to(self.device)
                
                next_state = self.model(state_history, current_history, state, current)
                
                # Shift history
                state_history[:, :-1] = state_history[:, 1:].clone()
                state_history[:, -1] = state
                
                current_history[:, :-1] = current_history[:, 1:].clone()
                current_history[:, -1] = current
                
                state = next_state
                trajectory.append(next_state.cpu().numpy()[0])
        
        return np.array(trajectory)
    
    def compute_rollout_error(self, test_trajectory: Dict,
                             steps_to_predict: int = 50) -> Dict:
        """Compute error on full trajectory prediction."""
        true_states = torch.from_numpy(test_trajectory['states']).float()
        true_currents = torch.from_numpy(test_trajectory['currents']).float()
        
        # Normalize using dataset stats
        true_states_norm = (true_states - torch.tensor(self.dataset.state_mean)) / \
                          torch.tensor(self.dataset.state_std)
        true_currents_norm = (true_currents - torch.tensor(self.dataset.currents_mean)) / \
                            torch.tensor(self.dataset.currents_std)
        
        initial_state = true_states_norm[0]
        current_sched = true_currents_norm[:steps_to_predict]
        
        pred_trajectory = self.predict_rollout(
            initial_state, current_sched, steps_to_predict
        )
        
        true_trajectory_norm = true_states_norm[:steps_to_predict].numpy()
        
        error = pred_trajectory - true_trajectory_norm
        
        return {
            'mae': np.mean(np.abs(error)),
            'rmse': np.sqrt(np.mean(error**2)),
            'max_error': np.max(np.abs(error)),
            'pred_trajectory': pred_trajectory,
            'true_trajectory': true_trajectory_norm,
            'error': error
        }
    
    def plot_rollout_prediction(self, test_trajectory: Dict, 
                               steps_to_predict: int = 50):
        """Plot predicted vs true trajectory."""
        metrics = self.compute_rollout_error(test_trajectory, steps_to_predict)
        
        pred = metrics['pred_trajectory']
        true = metrics['true_trajectory']
        
        fig, axes = plt.subplots(3, 1, figsize=(12, 8))
        
        for dim in range(3):
            axes[dim].plot(true[:, dim], 'b-', label='True', linewidth=2)
            axes[dim].plot(pred[:, dim], 'r--', label='Predicted', linewidth=2)
            axes[dim].set_ylabel(f'Dim {dim}')
            axes[dim].legend()
            axes[dim].grid()
        
        axes[-1].set_xlabel('Time step')
        plt.suptitle(f"Rollout Prediction (MAE={metrics['mae']:.4f}, "
                     f"RMSE={metrics['rmse']:.4f})")
        plt.tight_layout()
        plt.savefig('dynamics_rollout_prediction.png', dpi=150)
        print("Rollout prediction plot saved to dynamics_rollout_prediction.png")


# ============================================================================
# 6. MAIN TRAINING SCRIPT
# ============================================================================

def main():
    """Complete dynamics training pipeline."""
    
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}\n")
    
    # ---- STEP 1: Collect trajectory data ----
    print("=== Collecting Trajectory Data ===")
    crm_dyn = crm_cpp.CosseratRodDynamics("config.json")
    
    collector = DynamicsDataCollector(crm_dyn, save_dir="./data")
    trajectories = collector.collect_multiple_trajectories(
        num_trajectories=50,
        trajectory_length=100,
        dt=0.01,
        current_scale=0.05
    )
    collector.save_trajectories(trajectories, "dynamics_trajectories.pkl")
    
    # Or load existing:
    # trajectories = collector.load_trajectories("dynamics_trajectories.pkl")
    
    # ---- STEP 2: Create datasets ----
    print("\n=== Creating Datasets ===")
    full_dataset = SequencePredictionDataset(
        trajectories,
        seq_history_length=5,
        normalize=True
    )
    
    split_idx = int(0.8 * len(full_dataset))
    train_dataset = torch.utils.data.Subset(full_dataset, range(split_idx))
    val_dataset = torch.utils.data.Subset(full_dataset, range(split_idx, len(full_dataset)))
    
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
    
    print(f"Train size: {len(train_dataset)}, Val size: {len(val_dataset)}")
    
    # ---- STEP 3: Train different models ----
    print("\n=== Training Models ===")
    
    models_to_train = {
        'LSTM': LSTMDynamics(state_dim=3, action_dim=3, hidden_dim=64),
        'RNN': SimpleRNNDynamics(state_dim=3, action_dim=3, hidden_dim=64),
        'NeuralODE': NeuralODEDynamics(state_dim=3, action_dim=3, hidden_dim=64, dt=0.01)
    }
    
    best_model_name = None
    best_val_loss = float('inf')
    
    for model_name, model in models_to_train.items():
        print(f"\n--- Training {model_name} ---")
        trainer = DynamicsTrainer(model, device=device, learning_rate=1e-3)
        trainer.train(train_loader, val_loader, num_epochs=50, patience=10)
        trainer.plot_loss()
        
        final_val_loss = trainer.val_loss_history[-1]
        if final_val_loss < best_val_loss:
            best_val_loss = final_val_loss
            best_model_name = model_name
        
        torch.save(model.state_dict(), f"dynamics_{model_name.lower()}.pth")
    
    # ---- STEP 4: Evaluate best model ----
    print(f"\n=== Evaluation (Best: {best_model_name}) ===")
    
    best_model = models_to_train[best_model_name]
    best_model.load_state_dict(torch.load(f"dynamics_{best_model_name.lower()}.pth"))
    
    evaluator = DynamicsEvaluator(best_model, full_dataset, device=device)
    
    # Test on a few trajectories
    for i in range(3):
        test_traj = trajectories[i]
        metrics = evaluator.compute_rollout_error(test_traj, steps_to_predict=50)
        print(f"\nTest {i+1}: MAE={metrics['mae']:.6f}, "
              f"RMSE={metrics['rmse']:.6f}, Max={metrics['max_error']:.6f}")
    
    # Plot example rollout
    evaluator.plot_rollout_prediction(trajectories[0], steps_to_predict=50)
    
    print("\n=== Dynamics training complete ===")


if __name__ == "__main__":
    main()
