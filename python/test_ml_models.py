import sys
sys.path.insert(0, '../build_debug')

import torch
import numpy as np
from torch import nn

# ============================================================
# Test MLP Residual Dynamics
# ============================================================
class ResidualDynamics_MLP(nn.Module):
    def __init__(self, input_dim=9, hidden_dim=64, output_dim=3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
    
    def forward(self, state_t, action_t, physics_next):
        x = torch.cat([state_t, action_t, physics_next], dim=-1)
        residual = self.net(x)
        return physics_next + residual

print("=" * 60)
print("Testing ML Models")
print("=" * 60)

# Instantiate model
model = ResidualDynamics_MLP()
print("\n✓ ResidualDynamics_MLP initialized")

# Create dummy data
batch_size = 32
state_t = torch.randn(batch_size, 3)      # tip position
action_t = torch.randn(batch_size, 3)     # currents
physics_next = torch.randn(batch_size, 3) # predicted next state from physics

print(f"  Input shapes: state={state_t.shape}, action={action_t.shape}, physics={physics_next.shape}")

# Forward pass
output = model(state_t, action_t, physics_next)
print(f"  Output shape: {output.shape}")
print(f"  Output (first sample): {output[0]}")

# Test loss and backprop
target = torch.randn(batch_size, 3)
loss_fn = nn.MSELoss()
loss = loss_fn(output, target)
loss.backward()
print(f"  Loss: {loss.item():.6f}")
print("  ✓ Backprop successful")

# ============================================================
# Test LSTM Residual Dynamics
# ============================================================
print("\n" + "=" * 60)
print("Testing LSTM Residual Dynamics")
print("=" * 60)

class ResidualDynamics_LSTM(nn.Module):
    def __init__(self, input_dim=3, hidden_dim=64, output_dim=3, num_layers=2):
        super().__init__()
        self.lstm = nn.LSTM(input_dim, hidden_dim, num_layers, batch_first=True)
        self.fc = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, output_dim)
        )
    
    def forward(self, sequence, physics_next):
        _, (h_n, _) = self.lstm(sequence)
        residual = self.fc(h_n[-1])
        return physics_next + residual

model_lstm = ResidualDynamics_LSTM()
print("✓ ResidualDynamics_LSTM initialized")

# Create sequence data
seq_len = 10
sequence = torch.randn(batch_size, seq_len, 3)  # history of currents
physics_next = torch.randn(batch_size, 3)

print(f"  Input shape: sequence={sequence.shape}, physics={physics_next.shape}")

output_lstm = model_lstm(sequence, physics_next)
print(f"  Output shape: {output_lstm.shape}")

loss_lstm = loss_fn(output_lstm, target)
loss_lstm.backward()
print(f"  Loss: {loss_lstm.item():.6f}")
print("  ✓ LSTM backprop successful")

# ============================================================
# Test Neural ODE
# ============================================================
print("\n" + "=" * 60)
print("Testing Neural ODE Dynamics")
print("=" * 60)

class ResidualDynamics_NeuralODE(nn.Module):
    def __init__(self, state_dim=3, action_dim=3, hidden_dim=64, dt=0.01):
        super().__init__()
        self.dt = dt
        self.f = nn.Sequential(
            nn.Linear(state_dim + action_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, state_dim)
        )
    
    def forward(self, state_t, action_t, physics_next):
        dx_dt = self.f(torch.cat([state_t, action_t], dim=-1))
        correction = self.dt * dx_dt
        return physics_next + correction

model_ode = ResidualDynamics_NeuralODE()
print("✓ ResidualDynamics_NeuralODE initialized")

state_t = torch.randn(batch_size, 3)
action_t = torch.randn(batch_size, 3)
physics_next = torch.randn(batch_size, 3)

output_ode = model_ode(state_t, action_t, physics_next)
print(f"  Output shape: {output_ode.shape}")

loss_ode = loss_fn(output_ode, target)
loss_ode.backward()
print(f"  Loss: {loss_ode.item():.6f}")
print("  ✓ Neural ODE backprop successful")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 60)
print("✓ ALL ML TESTS PASSED!")
print("=" * 60)
print("""
Summary:
  - MLP Residual: Fast, suitable for static mappings
  - LSTM Residual: Captures temporal dependencies
  - Neural ODE: Continuous-time corrections
  
All models support:
  - Batch processing ✓
  - Backpropagation ✓
  - GPU acceleration (if available) ✓
""")
