#!/bin/bash
# Automated data generation

set -e

echo "=========================================="
echo "Generating Training Data"
echo "=========================================="

# Load config
if [ ! -f config.yaml ]; then
    echo "Error: config.yaml not found. Run ./scripts/setup.sh first"
    exit 1
fi

# Parse config (simple grep-based parsing)
CONFIG_DIM=$(grep "config_dim:" config.yaml | awk '{print $2}')
STATE_DIM=$(grep "state_dim:" config.yaml | awk '{print $2}')
N_SAMPLES=$(grep "n_samples:" config.yaml | awk '{print $2}')
N_TRAJ=$(grep "n_trajectories:" config.yaml | awk '{print $2}')
TRAJ_LEN=$(grep "trajectory_length:" config.yaml | awk '{print $2}')

echo "Configuration:"
echo "  Config dimension: $CONFIG_DIM"
echo "  State dimension: $STATE_DIM"
echo "  Kinematics samples: $N_SAMPLES"
echo "  Dynamics trajectories: $N_TRAJ"
echo ""

# Generate data
python3 examples/01_generate_data.py \
    --n_samples $N_SAMPLES \
    --n_trajectories $N_TRAJ \
    --trajectory_length $TRAJ_LEN \
    --output_dir data/ \
    --config_dim $CONFIG_DIM \
    --state_dim $STATE_DIM

echo ""
echo "✓ Data generation complete!"
echo "  Kinematics: data/kinematics_data.npz"
echo "  Dynamics: data/dynamics_data.npz"
