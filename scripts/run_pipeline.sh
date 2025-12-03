#!/bin/bash
# Complete automated pipeline

set -e

START_TIME=$(date +%s)

echo "=========================================="
echo "Running Complete ML/RL Pipeline"
echo "=========================================="
echo ""

# Step 1: Generate data
echo "Step 1/4: Generating data..."
./scripts/generate_data.sh

# Step 2: Train models
echo ""
echo "Step 2/4: Training models..."
./scripts/train_all.sh

# Step 3: Evaluate
echo ""
echo "Step 3/4: Evaluating agent..."
./scripts/evaluate.sh

# Step 4: Run tests
echo ""
echo "Step 4/4: Running tests..."
./scripts/test.sh

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "=========================================="
echo "✓ Pipeline Complete!"
echo "=========================================="
echo "Total time: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo ""
echo "Results:"
echo "  Data: data/"
echo "  Models: models/"
echo "  Evaluation: evaluation_trajectory.png"
echo "  Training curves: training_curves.png"
echo "  Test coverage: htmlcov/index.html"
