#!/bin/bash
# Complete automated setup script

set -e  # Exit on any error

echo "=========================================="
echo "MRI Catheter ML/RL Setup Script"
echo "=========================================="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
CRM_PATH="${CRM_PATH:-../CRM}"
CONFIG_DIM="${CONFIG_DIM:-6}"
STATE_DIM="${STATE_DIM:-24}"

echo -e "${YELLOW}Configuration:${NC}"
echo "  CRM_PATH: $CRM_PATH"
echo "  CONFIG_DIM: $CONFIG_DIM"
echo "  STATE_DIM: $STATE_DIM"
echo ""

# Function to print status
print_status() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_info() {
    echo -e "${YELLOW}→${NC} $1"
}

# Step 1: Check CRM code exists
echo "=========================================="
echo "Step 1: Checking CRM code"
echo "=========================================="

if [ ! -d "$CRM_PATH" ]; then
    print_error "CRM code not found at: $CRM_PATH"
    echo ""
    echo "Please set CRM_PATH environment variable:"
    echo "  export CRM_PATH=/path/to/your/CRM"
    echo "  ./scripts/setup.sh"
    exit 1
fi

if [ ! -f "$CRM_PATH/CMakeLists.txt" ]; then
    print_error "CRM CMakeLists.txt not found"
    exit 1
fi

print_status "CRM code found at: $CRM_PATH"
echo ""

# Step 2: Build CRM libraries if needed
echo "=========================================="
echo "Step 2: Building CRM libraries"
echo "=========================================="

CRM_BUILD_DIR="$CRM_PATH/build"

if [ ! -f "$CRM_BUILD_DIR/libCRMCPPLib.a" ]; then
    print_info "CRM libraries not found, building..."
    
    mkdir -p "$CRM_BUILD_DIR"
    cd "$CRM_BUILD_DIR"
    
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j$(nproc)
    
    cd - > /dev/null
    
    if [ ! -f "$CRM_BUILD_DIR/libCRMCPPLib.a" ]; then
        print_error "Failed to build libCRMCPPLib.a"
        exit 1
    fi

    print_status "CRM libraries built successfully"
else
    print_status "CRM libraries already built"
fi
echo ""

# Step 3: Build C++ bindings
echo "=========================================="
echo "Step 3: Building C++ Python bindings"
echo "=========================================="

BINDINGS_DIR="crm_ml_rl/bindings"
BINDINGS_BUILD_DIR="$BINDINGS_DIR/build"

mkdir -p "$BINDINGS_BUILD_DIR"
cd "$BINDINGS_BUILD_DIR"

print_info "Running CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCRM_INCLUDE_DIR="$(realpath $CRM_PATH)/src" \
    -DCRM_LIB_DIR="$(realpath $CRM_PATH)/build"

print_info "Building bindings (using $(nproc) cores)..."
make -j$(nproc)

cd - > /dev/null

# Check if binding was created
BINDING_FILE=$(find crm_ml_rl -name "crm_cpp*.so" 2>/dev/null | head -1)

if [ -z "$BINDING_FILE" ]; then
    print_error "Failed to build C++ bindings"
    exit 1
fi

print_status "C++ bindings built: $BINDING_FILE"
echo ""

# Step 4: Install Python package
echo "=========================================="
echo "Step 4: Installing Python package"
echo "=========================================="

print_info "Installing in editable mode..."
pip install -e . > /dev/null 2>&1

# Verify installation
if python3 -c "import crm_cpp" 2>/dev/null; then
    print_status "C++ bindings importable"
else
    print_error "Failed to import C++ bindings"
    exit 1
fi

if python3 -c "from crm_ml_rl import ForwardKinematicsAPI" 2>/dev/null; then
    print_status "Python package installed successfully"
else
    print_error "Failed to import Python package"
    exit 1
fi
echo ""

# Step 4.5: Check GPU and offer to enable
echo "=========================================="
echo "Step 4.5: Checking GPU support"
echo "=========================================="

if python3 -c "import torch; print(torch.cuda.is_available())" 2>/dev/null | grep -q "True"; then
    print_status "GPU support already enabled"
elif command -v nvidia-smi &> /dev/null; then
    print_info "NVIDIA GPU detected but PyTorch is CPU-only"
    echo ""
    read -p "Enable GPU support? (Y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Nn]$ ]]; then
        ./scripts/enable_gpu.sh
    fi
else
    print_info "No GPU detected, using CPU (slower but works)"
fi
echo ""

# Step 5: Create necessary directories
echo "=========================================="
echo "Step 5: Creating directories"
echo "=========================================="

mkdir -p data
mkdir -p models
mkdir -p logs
mkdir -p results

print_status "Created data/, models/, logs/, results/"
echo ""

# Step 6: Save configuration
echo "=========================================="
echo "Step 6: Saving configuration"
echo "=========================================="

cat > config.yaml <<EOF
# Auto-generated configuration
catheter:
  config_dim: $CONFIG_DIM
  state_dim: $STATE_DIM
  
paths:
  crm_path: $CRM_PATH
  data_dir: data
  models_dir: models
  logs_dir: logs
  results_dir: results

training:
  n_samples: 50000
  n_trajectories: 1000
  trajectory_length: 50
  
  kinematics:
    epochs: 100
    batch_size: 512
    lr: 0.001
  
  dynamics:
    epochs: 100
    batch_size: 256
    lr: 0.001
  
  rl:
    n_episodes: 5000
    batch_size: 256
    lr: 0.0003
    
evaluation:
  n_trials: 100
EOF

print_status "Configuration saved to config.yaml"
echo ""

# Done!
echo "=========================================="
echo -e "${GREEN}Setup Complete!${NC}"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Generate training data:"
echo "     ./scripts/generate_data.sh"
echo ""
echo "  2. Train all models:"
echo "     ./scripts/train_all.sh"
echo ""
echo "  3. Run full pipeline:"
echo "     ./scripts/run_pipeline.sh"
echo ""
echo "  4. Run tests:"
echo "     ./scripts/test.sh"
echo ""
