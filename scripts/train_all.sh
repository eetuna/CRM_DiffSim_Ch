#!/bin/bash
# Train all models sequentially

set -e

echo "=========================================="
echo "Training All Models"
echo "=========================================="

# Load config
if [ ! -f config.yaml ]; then
    echo "Error: config.yaml not found. Run ./scripts/setup.sh first"
    exit 1
fi

# Parse training config
KIN_EPOCHS=$(grep -A 3 "kinematics:" config.yaml | grep "epochs:" | awk '{print $2}')
KIN_BATCH=$(grep -A 3 "kinematics:" config.yaml | grep "batch_size:" | awk '{print $2}')
KIN_LR=$(grep -A 3 "kinematics:" config.yaml | grep "lr:" | awk '{print $2}')

DYN_EPOCHS=$(grep -A 3 "dynamics:" config.yaml | grep "epochs:" | awk '{print $2}')
DYN_BATCH=$(grep -A 3 "dynamics:" config.yaml | grep "batch_size:" | awk '{print $2}')
DYN_LR=$(grep -A 3 "dynamics:" config.yaml | grep "lr:" | awk '{print $2}')

RL_EPISODES=$(grep -A 4 "rl:" config.yaml | grep "n_episodes:" | awk '{print $2}')

echo ""
echo "=========================================="
echo "Step 1/3: Training Kinematics Model"
echo "=========================================="
python3 examples/02_train_kinematics.py \
    --data_path data/kinematics_data.npz \
    --epochs $KIN_EPOCHS \
    --batch_size $KIN_BATCH \
    --lr $KIN_LR \
    --save_path models/kinematics_residual.pth

echo ""
echo "=========================================="
echo "Step 2/3: Training Dynamics Model"
echo "=========================================="
python3 examples/03_train_dynamics.py \
    --data_path data/dynamics_data.npz \
    --epochs $DYN_EPOCHS \
    --batch_size $DYN_BATCH \
    --lr $DYN_LR \
    --save_path models/dynamics_model.pth

echo ""
echo "=========================================="
echo "Step 3/3: Training RL Agent"
echo "=========================================="
python3 examples/04_train_rl.py \
    --n_episodes $RL_EPISODES \
    --dynamics_model models/dynamics_model.pth \
    --save_path models/td3_agent.pth

echo ""
echo "=========================================="
echo "✓ All models trained successfully!"
echo "=========================================="
echo "  Kinematics: models/kinematics_residual.pth"
echo "  Dynamics: models/dynamics_model.pth"
echo "  RL Agent: models/td3_agent.pth"
