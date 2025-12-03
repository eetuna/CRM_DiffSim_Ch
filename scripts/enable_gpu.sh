#!/bin/bash
# Enable GPU support by installing CUDA-enabled PyTorch

set -e

echo "=========================================="
echo "GPU Support Configuration"
echo "=========================================="
echo ""

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_status() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_info() {
    echo -e "${YELLOW}→${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}⚠${NC} $1"
}

# Check current PyTorch installation
print_info "Checking current PyTorch installation..."
CURRENT_TORCH=$(python3 -c "import torch; print(torch.__version__)" 2>/dev/null || echo "not installed")
CUDA_AVAILABLE=$(python3 -c "import torch; print(torch.cuda.is_available())" 2>/dev/null || echo "False")

echo "  Current PyTorch: $CURRENT_TORCH"
echo "  CUDA available: $CUDA_AVAILABLE"
echo ""

# Check if NVIDIA GPU exists
print_info "Checking for NVIDIA GPU..."
if command -v nvidia-smi &> /dev/null; then
    GPU_INFO=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "Detection failed")
    print_status "NVIDIA GPU detected: $GPU_INFO"
    
    CUDA_VERSION=$(nvidia-smi | grep "CUDA Version" | awk '{print $9}' || echo "unknown")
    echo "  CUDA Driver Version: $CUDA_VERSION"
    echo ""
else
    print_warning "No NVIDIA GPU detected (nvidia-smi not found)"
    echo ""
    echo "This script will install CUDA PyTorch, but it won't help without a GPU."
    echo "Training will fall back to CPU."
    echo ""
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
fi

# Determine CUDA version to install
if [[ "$CUDA_VERSION" == "unknown" ]] || [[ "$CUDA_VERSION" == "" ]]; then
    CUDA_INSTALL="cu118"  # Default to CUDA 11.8
    print_warning "Could not detect CUDA version, defaulting to CUDA 11.8"
else
    CUDA_MAJOR=$(echo $CUDA_VERSION | cut -d. -f1)
    CUDA_MINOR=$(echo $CUDA_VERSION | cut -d. -f2)
    
    if [[ "$CUDA_MAJOR" == "12" ]]; then
        CUDA_INSTALL="cu121"
        print_info "Detected CUDA 12.x, will install PyTorch with CUDA 12.1"
    elif [[ "$CUDA_MAJOR" == "11" ]]; then
        CUDA_INSTALL="cu118"
        print_info "Detected CUDA 11.x, will install PyTorch with CUDA 11.8"
    else
        CUDA_INSTALL="cu118"
        print_warning "Unsupported CUDA version, defaulting to CUDA 11.8"
    fi
fi

echo ""
echo "=========================================="
echo "Installing GPU-Enabled PyTorch"
echo "=========================================="
echo ""

# Uninstall CPU-only PyTorch
print_info "Uninstalling CPU-only PyTorch..."
pip uninstall -y torch torchvision torchaudio 2>/dev/null || true

# Install GPU-enabled PyTorch
print_info "Installing PyTorch with CUDA $CUDA_INSTALL..."
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/$CUDA_INSTALL

# Verify installation
echo ""
print_info "Verifying GPU support..."
python3 << 'EOF'
import torch
print(f"PyTorch version: {torch.__version__}")
print(f"CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"CUDA version: {torch.version.cuda}")
    print(f"GPU count: {torch.cuda.device_count()}")
    for i in range(torch.cuda.device_count()):
        print(f"  GPU {i}: {torch.cuda.get_device_name(i)}")
    print(f"Current device: {torch.cuda.current_device()}")
else:
    print("WARNING: CUDA not available!")
EOF

# Check status
if python3 -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
    echo ""
    echo "=========================================="
    echo -e "${GREEN}GPU Support Enabled!${NC}"
    echo "=========================================="
    echo ""
    
    # Show memory info
    print_info "GPU Memory Info:"
    nvidia-smi --query-gpu=memory.total,memory.free,memory.used --format=csv
    
    echo ""
    echo "Your training will now use GPU acceleration!"
    echo ""
    echo "Tips:"
    echo "  - Monitor GPU usage: watch -n 1 nvidia-smi"
    echo "  - Training will be 5-20x faster than CPU"
    echo "  - Increase batch size to utilize GPU better"
    echo ""
else
    echo ""
    echo "=========================================="
    echo -e "${YELLOW}Warning: GPU Support Not Available${NC}"
    echo "=========================================="
    echo ""
    echo "PyTorch installed but CUDA not available."
    echo "Training will use CPU (slower)."
    echo ""
    echo "Possible reasons:"
    echo "  - No NVIDIA GPU in system"
    echo "  - CUDA drivers not installed"
    echo "  - Driver/PyTorch version mismatch"
    echo ""
fi
