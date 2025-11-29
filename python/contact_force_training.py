import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
import numpy as np
import crm_cpp

class ContactDataCollector:
    """Collect contact force data from your contact force Jacobian."""
    
    def __init__(self, cfj_model: crm_cpp.ContactForceJacobian):
        self.cfj = cfj_model
        self.data = {'currents': [], 'jacobian': [], 'contact_force': []}
    
    def collect_contact_data(self, num_samples: int = 200, 
                             contact_surface_samples: int = 5):
        """Collect contact scenarios."""
        print(f"Collecting {num_samples} contact samples...")
        
        for i in range(num_samples):
            currents = np.random.uniform(-0.1, 0.1, 3)
            inserted_length = np.random.uniform(0.07, 0.10)
            
            # Sample contact positions on known surface
            for _ in range(contact_surface_samples):
                contact_pos = np.random.randn(3) * 0.01  # Small perturbations
                
                jacobian = np.array(
                    self.cfj.compute_jacobian(
                        currents.tolist(), inserted_length, contact_pos.tolist()
                    )
                )
                
                # Simulate contact force (in practice, measured)
                delta_currents = np.random.randn(3) * 0.01
                contact_force = np.random.randn(3) * 0.001  # Small forces
                
                self.data['currents'].append(currents)
                self.data['jacobian'].append(jacobian)
                self.data['contact_force'].append(contact_force)
        
        for key in self.data:
            self.data[key] = np.array(self.data[key])
    
    def save(self, filename: str = "contact_data.pkl"):
        import pickle
        with open(filename, 'wb') as f:
            pickle.dump(self.data, f)


class ContactForceDataset(Dataset):
    """Dataset: [currents + jacobian_flattened] → contact_force."""
    
    def __init__(self, currents: np.ndarray, jacobians: np.ndarray, 
                 forces: np.ndarray):
        # Flatten jacobians: (N, 3, 3) → (N, 9)
        jacobians_flat = jacobians.reshape(jacobians.shape[0], -1)
        self.x = torch.from_numpy(
            np.concatenate([currents, jacobians_flat], axis=1)
        ).float()
        self.y = torch.from_numpy(forces).float()
    
    def __len__(self):
        return len(self.x)
    
    def __getitem__(self, idx):
        return {'x': self.x[idx], 'y': self.y[idx]}


def train_contact_force_model():
    """Train contact force predictor."""
    
    cfj = crm_cpp.ContactForceJacobian("config.json")
    collector = ContactDataCollector(cfj)
    collector.collect_contact_data(num_samples=200)
    
    # Create dataset
    split_idx = int(0.8 * len(collector.data['currents']))
    
    train_data = ContactForceDataset(
        collector.data['currents'][:split_idx],
        collector.data['jacobian'][:split_idx],
        collector.data['contact_force'][:split_idx]
    )
    val_data = ContactForceDataset(
        collector.data['currents'][split_idx:],
        collector.data['jacobian'][split_idx:],
        collector.data['contact_force'][split_idx:]
    )
    
    train_loader = DataLoader(train_data, batch_size=16, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=16, shuffle=False)
    
    # Train model
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = nn.Sequential(
        nn.Linear(12, 128),
        nn.ReLU(),
        nn.Dropout(0.2),
        nn.Linear(128, 128),
        nn.ReLU(),
        nn.Dropout(0.2),
        nn.Linear(128, 3)
    ).to(device)
    
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.MSELoss()
    
    for epoch in range(30):
        model.train()
        for batch in train_loader:
            x = batch['x'].to(device)
            y = batch['y'].to(device)
            pred = model(x)
            loss = criterion(pred, y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        
        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch in val_loader:
                x = batch['x'].to(device)
                y = batch['y'].to(device)
                pred = model(x)
                val_loss += criterion(pred, y).item()
        
        print(f"Epoch {epoch+1}: Train Loss={loss.item():.6f}, "
              f"Val Loss={val_loss/len(val_loader):.6f}")
    
    torch.save(model.state_dict(), "contact_force_model.pth")
    print("Contact force model saved.")


if __name__ == "__main__":
    train_contact_force_model()
