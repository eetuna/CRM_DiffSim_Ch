import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import matplotlib.pyplot as plt
from pathlib import Path
import pickle
import json
from typing import Tuple, List, Dict, Optional
import sys

sys.path.insert(0, '../build_debug')
import crm_cpp  # Your C++ binding


# ============================================================================
# 1. DATA COLLECTION & PREPROCESSING
# ============================================================================

class CatheterDataCollector:
    """Collects experimental or simulated data from camera + physics model."""

    def __init__(self, config_file: str, save_dir: str = "../data"):
        """
        Parameters
        ----------
        config_file : str
            Path to catheter configuration JSON used by CRM_ForwardKinematics.
        save_dir : str
            Directory to store collected datasets.
        """
        self.config_file = config_file
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(exist_ok=True)
        self.data = {
            'currents': [],
            'inserted_length': [],
            'model_tip': [],
            'camera_tip': [],
            'residual': []
        }

    def _fk(self, currents: np.ndarray, inserted_length: float) -> np.ndarray:
        """
        Call C++ forward kinematics and return only tip position (first 3 entries).
        """
        pRu = np.array(
            crm_cpp.forward_kinematics(currents, inserted_length, self.config_file)
        )
        return pRu[:3]  # p_x, p_y, p_z

    def collect_free_space_data(self, num_samples: int = 500,
                                camera_noise_std: float = 0.5e-3):
        """
        Sweep currents, collect (currents, model_pred, camera_measured) triplets.

        In practice, replace camera_measured with real camera data.
        """
        print(f"Collecting {num_samples} free-space samples...")

        # Define current sweep ranges (adjust to your system)
        n_coils = crm_cpp.get_num_actuators()
        current_ranges = [(-0.1, 0.1)] * n_coils
        inserted_length_range = (0.07, 0.10)  # 70-100 mm, example

        for i in range(num_samples):
            currents = np.array([
                np.random.uniform(lo, hi)
                for lo, hi in current_ranges
            ], dtype=np.float64)
            inserted_length = float(np.random.uniform(*inserted_length_range))

            # Physics model prediction
            model_tip = self._fk(currents, inserted_length)

            # Simulated camera measurement
            camera_tip = model_tip + np.random.randn(3) * camera_noise_std

            residual = camera_tip - model_tip

            self.data['currents'].append(currents)
            self.data['inserted_length'].append(inserted_length)
            self.data['model_tip'].append(model_tip)
            self.data['camera_tip'].append(camera_tip)
            self.data['residual'].append(residual)

            if (i + 1) % 100 == 0:
                print(f"  {i + 1}/{num_samples} samples collected")

        # Convert to numpy
        for key in self.data:
            self.data[key] = np.array(self.data[key])

        print(f"Data shape: currents {self.data['currents'].shape}, "
              f"residuals {self.data['residual'].shape}")

    def normalize_data(self) -> Dict:
        """Compute and store normalization statistics."""
        norm_stats = {
            'currents_mean': self.data['currents'].mean(axis=0),
            'currents_std': self.data['currents'].std(axis=0),
            'residual_mean': self.data['residual'].mean(axis=0),
            'residual_std': self.data['residual'].std(axis=0),
        }

        # Normalize
        self.data['currents_norm'] = (
            (self.data['currents'] - norm_stats['currents_mean']) /
            (norm_stats['currents_std'] + 1e-8)
        )
        self.data['residual_norm'] = (
            (self.data['residual'] - norm_stats['residual_mean']) /
            (norm_stats['residual_std'] + 1e-8)
        )

        return norm_stats

    def save_data(self, filename: str = "catheter_data.pkl"):
        """Save collected data."""
        filepath = self.save_dir / filename
        with open(filepath, 'wb') as f:
            pickle.dump(self.data, f)
        print(f"Data saved to {filepath}")

    def load_data(self, filename: str = "catheter_data.pkl"):
        """Load previously collected data."""
        filepath = self.save_dir / filename
        with open(filepath, 'rb') as f:
            self.data = pickle.load(f)
        print(f"Data loaded from {filepath}")


# ============================================================================
# 2. DATASET CLASSES FOR PyTorch
# ============================================================================

class ResidualDataset(Dataset):
    """Dataset for residual learning: (currents, model_tip) → residual."""

    def __init__(self, currents: np.ndarray, model_tip: np.ndarray,
                 residual: np.ndarray, norm_stats: Optional[Dict] = None):
        self.currents = torch.from_numpy(currents).float()
        self.model_tip = torch.from_numpy(model_tip).float()
        self.residual = torch.from_numpy(residual).float()
        self.norm_stats = norm_stats

        # Normalize currents if stats provided
        if norm_stats:
            self.currents = (
                (self.currents - torch.tensor(norm_stats['currents_mean']).float()) /
                (torch.tensor(norm_stats['currents_std']).float() + 1e-8)
            )

    def __len__(self):
        return len(self.currents)

    def __getitem__(self, idx):
        return {
            'currents': self.currents[idx],
            'model_tip': self.model_tip[idx],
            'residual': self.residual[idx]
        }


class DynamicsDataset(Dataset):
    """Dataset for dynamics: (currents_t, state_t) → state_t+1."""

    def __init__(self, sequences: List[np.ndarray], seq_length: int = 10):
        """
        sequences: list of (N, 6) trajectories [currents (3) + tip_pos (3)]
        """
        self.data = []
        for seq in sequences:
            for i in range(len(seq) - seq_length):
                self.data.append({
                    'history': seq[i:i + seq_length],   # (seq_length, 6)
                    'next_state': seq[i + seq_length]   # (6,)
                })

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return {
            'history': torch.from_numpy(self.data[idx]['history']).float(),
            'next_state': torch.from_numpy(self.data[idx]['next_state']).float()
        }


# ============================================================================
# 3. NEURAL NETWORK ARCHITECTURES
# ============================================================================

class ResidualMLP(nn.Module):
    """Small MLP to learn residuals: input=[currents, model_tip] → output=residual."""

    def __init__(self, input_dim: int = 6, hidden_dim: int = 64, output_dim: int = 3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x):
        return self.net(x)


class DynamicsLSTM(nn.Module):
    """LSTM for learning dynamics (if needed for history-dependent effects)."""

    def __init__(self, input_dim: int = 6, hidden_dim: int = 64, output_dim: int = 6):
        super().__init__()
        self.lstm = nn.LSTM(input_size=input_dim, hidden_size=hidden_dim,
                            num_layers=2, batch_first=True)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x):
        """x: (batch, seq_length, input_dim)"""
        _, (h_n, _) = self.lstm(x)
        out = self.fc(h_n[-1])  # last layer hidden state
        return out


class ContactForcePredictor(nn.Module):
    """Learn contact force residuals on top of analytical Jacobian."""

    def __init__(self, input_dim: int = 12, hidden_dim: int = 128, output_dim: int = 3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x):
        return self.net(x)


# ============================================================================
# 4. TRAINING LOOP
# ============================================================================

class ModelTrainer:
    """Generic trainer for any PyTorch model."""

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
        self.model.train()
        total_loss = 0.0

        for batch in train_loader:
            if isinstance(batch, dict):
                batch = {k: v.to(self.device) for k, v in batch.items()}

            if 'history' in batch:  # Dynamics model
                pred = self.model(batch['history'])
                target = batch['next_state']
            else:  # Residual model
                x = torch.cat([batch['currents'], batch['model_tip']], dim=-1)
                pred = self.model(x)
                target = batch['residual']

            loss = self.criterion(pred, target)

            self.optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.optimizer.step()

            total_loss += loss.item()

        avg_loss = total_loss / len(train_loader)
        self.train_loss_history.append(avg_loss)
        return avg_loss

    def validate(self, val_loader: DataLoader) -> float:
        self.model.eval()
        total_loss = 0.0

        with torch.no_grad():
            for batch in val_loader:
                if isinstance(batch, dict):
                    batch = {k: v.to(self.device) for k, v in batch.items()}

                if 'history' in batch:
                    pred = self.model(batch['history'])
                    target = batch['next_state']
                else:
                    x = torch.cat([batch['currents'], batch['model_tip']], dim=-1)
                    pred = self.model(x)
                    target = batch['residual']

                loss = self.criterion(pred, target)
                total_loss += loss.item()

        avg_loss = total_loss / len(val_loader)
        self.val_loss_history.append(avg_loss)
        return avg_loss

    def train(self, train_loader: DataLoader, val_loader: DataLoader,
              num_epochs: int = 50, patience: int = 10):
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
                self.save_checkpoint("best_model.pth")
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
        plt.savefig('loss_curve.png', dpi=150)
        print("Loss curve saved to loss_curve.png")


# ============================================================================
# 5. INFERENCE & EVALUATION
# ============================================================================

class CatheterPredictor:
    """Unified predictor: uses physics model + learned residuals."""

    def __init__(self, config_file: str,
                 residual_net: nn.Module, norm_stats: Dict,
                 device: str = 'cpu'):
        self.config_file = config_file
        self.residual_net = residual_net.to(device).eval()
        self.norm_stats = norm_stats
        self.device = device

    def _fk(self, currents: np.ndarray, inserted_length: float) -> np.ndarray:
        pRu = np.array(
            crm_cpp.forward_kinematics(currents, inserted_length, self.config_file)
        )
        return pRu[:3]

    def predict_tip(self, currents: np.ndarray, inserted_length: float) -> np.ndarray:
        """
        Predict tip position:
        tip_final = crm_prediction + residual_net_correction
        """
        model_tip = self._fk(currents, inserted_length)

        currents_norm = (
            (currents - self.norm_stats['currents_mean']) /
            (self.norm_stats['currents_std'] + 1e-8)
        )

        x = np.concatenate([currents_norm, model_tip])
        x_tensor = torch.from_numpy(x).float().to(self.device).unsqueeze(0)

        with torch.no_grad():
            residual_pred = self.residual_net(x_tensor).cpu().numpy()[0]

        residual_pred = (
            residual_pred * (self.norm_stats['residual_std'] + 1e-8) +
            self.norm_stats['residual_mean']
        )

        return model_tip + residual_pred

    def compute_error_metrics(self, currents: np.ndarray,
                              inserted_length: float,
                              true_tip: np.ndarray) -> Dict:
        pred_tip = self.predict_tip(currents, inserted_length)
        error = pred_tip - true_tip

        return {
            'mae': float(np.mean(np.abs(error))),
            'rmse': float(np.sqrt(np.mean(error ** 2))),
            'max_error': float(np.max(np.abs(error))),
            'error_vector': error,
        }


# ============================================================================
# 6. MAIN TRAINING SCRIPT
# ============================================================================

def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")

    # ---- STEP 1: Config file path ----
    print("\n=== Using CRM config ===")
    config_file = "../catheterdata/CatheterSpatialConfiguration_1.txt"  # adjust path as needed

    # ---- STEP 2: Collect or load data ----
    print("\n=== Data Collection ===")
    collector = CatheterDataCollector(config_file, save_dir="../data")

    collector.collect_free_space_data(num_samples=500, camera_noise_std=0.5e-3)
    norm_stats = collector.normalize_data()
    collector.save_data("catheter_data.pkl")

    # ---- STEP 3: Create datasets and dataloaders ----
    print("\n=== Creating Datasets ===")

    N = len(collector.data['currents'])
    split_idx = int(0.8 * N)

    train_dataset = ResidualDataset(
        collector.data['currents'][:split_idx],
        collector.data['model_tip'][:split_idx],
        collector.data['residual'][:split_idx],
        norm_stats=norm_stats,
    )

    val_dataset = ResidualDataset(
        collector.data['currents'][split_idx:],
        collector.data['model_tip'][split_idx:],
        collector.data['residual'][split_idx:],
        norm_stats=norm_stats,
    )

    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)

    # ---- STEP 4: Train residual model ----
    print("\n=== Training Residual MLP ===")

    # Input dim = currents_norm (n_coils) + model_tip (3)
    n_coils = crm_cpp.get_num_actuators()
    input_dim = n_coils + 3
    residual_model = ResidualMLP(input_dim=input_dim, hidden_dim=64, output_dim=3)
    trainer = ModelTrainer(residual_model, device=device, learning_rate=1e-3)
    trainer.train(train_loader, val_loader, num_epochs=50, patience=10)
    trainer.plot_loss()

    # ---- STEP 5: Evaluation ----
    print("\n=== Evaluation ===")
    trainer.load_checkpoint("best_model.pth")

    predictor = CatheterPredictor(config_file, residual_model, norm_stats, device=device)

    for i in range(5):
        idx = split_idx + i
        currents = collector.data['currents'][idx]
        true_tip = collector.data['camera_tip'][idx]
        inserted_length = collector.data['inserted_length'][idx]

        metrics = predictor.compute_error_metrics(currents, inserted_length, true_tip)
        print(f"Test {i+1}: "
              f"MAE={metrics['mae']*1e3:.3f} mm, "
              f"RMSE={metrics['rmse']*1e3:.3f} mm, "
              f"Max={metrics['max_error']*1e3:.3f} mm")

    torch.save(residual_model.state_dict(), "residual_model.pth")
    with open("norm_stats.json", 'w') as f:
        json.dump({k: v.tolist() for k, v in norm_stats.items()}, f)
    print("\nModel and normalization stats saved.")


if __name__ == "__main__":
    main()
