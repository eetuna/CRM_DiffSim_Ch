#!/bin/bash
# Automated evaluation

set -e

echo "=========================================="
echo "Evaluating RL Agent"
echo "=========================================="

# Check if agent exists
if [ ! -f models/td3_agent.pth ]; then
    echo "Error: Agent not found. Train first with ./scripts/train_all.sh"
    exit 1
fi

# Load config
N_TRIALS=$(grep "n_trials:" config.yaml | awk '{print $2}')
STATE_DIM=$(grep "state_dim:" config.yaml | awk '{print $2}')

python3 examples/05_evaluate.py \
    --agent models/td3_agent.pth \
    --dynamics_model models/dynamics_model.pth \
    --n_trials $N_TRIALS \
    --state_dim $STATE_DIM

echo ""
echo "✓ Evaluation complete!"
echo "  Results: results/"
echo "  Trajectory: evaluation_trajectory.png"
