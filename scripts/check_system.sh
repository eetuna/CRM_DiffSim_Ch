#!/bin/bash
# Check system requirements and configuration

echo "=========================================="
echo "System Check"
echo "=========================================="
echo ""

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

check_ok() {
    echo -e "${GREEN}✓${NC} $1"
}

check_fail() {
    echo -e "${RED}✗${NC} $1"
}

check_warn() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Python
echo "Python:"
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version | awk '{print $2}')
    check_ok "Python $PYTHON_VERSION"
else
    check_fail "Python 3 not found"
fi

# CMake
echo ""
echo "Build tools:"
if command -v cmake &> /dev/null; then
    CMAKE_VERSION=$(cmake --version | head -1 | awk '{print $3}')
    check_ok "CMake $CMAKE_VERSION"
else
    check_fail "CMake not found"
fi

# C++ compiler
if command -v g++ &> /dev/null; then
    GCC_VERSION=$(g++ --version | head -1 | awk '{print $3}')
    check_ok "g++ $GCC_VERSION"
else
    check_fail "g++ not found"
fi

# Eigen3
if pkg-config --exists eigen3; then
    EIGEN_VERSION=$(pkg-config --modversion eigen3)
    check_ok "Eigen3 $EIGEN_VERSION"
else
    check_warn "Eigen3 not found via pkg-config (may still be installed)"
fi

# PyTorch
echo ""
echo "PyTorch:"
if python3 -c "import torch" 2>/dev/null; then
    TORCH_VERSION=$(python3 -c "import torch; print(torch.__version__)")
    CUDA_AVAILABLE=$(python3 -c "import torch; print(torch.cuda.is_available())")
    check_ok "PyTorch $TORCH_VERSION"
    
    if [ "$CUDA_AVAILABLE" == "True" ]; then
        CUDA_VERSION=$(python3 -c "import torch; print(torch.version.cuda)")
        GPU_COUNT=$(python3 -c "import torch; print(torch.cuda.device_count())")
        check_ok "CUDA $CUDA_VERSION ($GPU_COUNT GPU(s))"
    else
        check_warn "CUDA not available (CPU-only)"
    fi
else
    check_fail "PyTorch not installed"
fi

# C++ bindings
echo ""
echo "C++ bindings:"
if python3 -c "import crm_cpp" 2>/dev/null; then
    check_ok "crm_cpp module built and importable"
else
    check_warn "crm_cpp module not found (run: make setup)"
fi

# Package installation
echo ""
echo "Python package:"
if python3 -c "from crm_ml_rl import ForwardKinematicsAPI" 2>/dev/null; then
    check_ok "crm_ml_rl package installed"
else
    check_warn "crm_ml_rl package not installed (run: pip install -e .)"
fi

# CRM code
echo ""
echo "CRM C++ code:"
CRM_PATH="${CRM_PATH:-../CRM}"
if [ -d "$CRM_PATH" ]; then
    check_ok "CRM path: $CRM_PATH"
    
    if [ -f "$CRM_PATH/build/libcrm.so" ]; then
        check_ok "libcrm.so found"
    else
        check_warn "libcrm.so not found (needs building)"
    fi
    
    if [ -f "$CRM_PATH/build/libcrmdyn.so" ]; then
        check_ok "libcrmdyn.so found"
    else
        check_warn "libcrmdyn.so not found (needs building)"
    fi
else
    check_fail "CRM code not found at: $CRM_PATH"
    echo "       Set CRM_PATH environment variable"
fi

# Directories
echo ""
echo "Project structure:"
for dir in data models logs results; do
    if [ -d "$dir" ]; then
        check_ok "$dir/ exists"
    else
        check_warn "$dir/ not found (will be created)"
    fi
done

# Config file
echo ""
echo "Configuration:"
if [ -f "config.yaml" ]; then
    check_ok "config.yaml exists"
else
    check_warn "config.yaml not found (run: make setup)"
fi

# GPU details (if available)
if command -v nvidia-smi &> /dev/null; then
    echo ""
    echo "GPU Information:"
    nvidia-smi --query-gpu=index,name,driver_version,memory.total --format=csv,noheader | while read line; do
        check_ok "GPU $line"
    done
fi

echo ""
echo "=========================================="
echo "System check complete"
echo "=========================================="
echo ""

# Recommendations
if ! python3 -c "import crm_cpp" 2>/dev/null; then
    echo "Next steps:"
    echo "  1. Run: make setup"
    echo "  2. Run: make pipeline"
fi

if python3 -c "import torch; assert not torch.cuda.is_available()" 2>/dev/null; then
    echo ""
    echo "Optional: Enable GPU support"
    echo "  Run: ./scripts/enable_gpu.sh"
fi
