#!/bin/bash
# Build C++ bindings only (useful for development)

set -e

echo "=========================================="
echo "Building C++ Bindings"
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

# Get CRM path
CRM_PATH="${CRM_PATH:-../CRM}"

if [ ! -d "$CRM_PATH" ]; then
    print_error "CRM code not found at: $CRM_PATH"
    echo ""
    echo "Set CRM_PATH environment variable:"
    echo "  export CRM_PATH=/path/to/your/CRM"
    echo "  ./scripts/build_cpp.sh"
    exit 1
fi

print_status "CRM path: $CRM_PATH"

# Check if CRM libraries exist
CRM_BUILD_DIR="$CRM_PATH/build"

if [ ! -f "$CRM_BUILD_DIR/libcrm.so" ]; then
    print_error "CRM libraries not found at: $CRM_BUILD_DIR"
    echo ""
    echo "Build your CRM libraries first:"
    echo "  cd $CRM_PATH"
    echo "  mkdir -p build && cd build"
    echo "  cmake .."
    echo "  make -j$(nproc)"
    exit 1
fi

print_status "CRM libraries found"
echo ""

# Clean previous build (optional)
if [ "$1" == "--clean" ]; then
    print_info "Cleaning previous build..."
    rm -rf crm_ml_rl/bindings/build/
    find crm_ml_rl -name "crm_cpp*.so" -delete
    print_status "Clean complete"
fi

# Build bindings
BINDINGS_DIR="crm_ml_rl/bindings"
BINDINGS_BUILD_DIR="$BINDINGS_DIR/build"

mkdir -p "$BINDINGS_BUILD_DIR"
cd "$BINDINGS_BUILD_DIR"

print_info "Configuring CMake..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCRM_INCLUDE_DIR="$(realpath $CRM_PATH)/include" \
    -DCRM_LIB_DIR="$(realpath $CRM_PATH)/build"

print_info "Building (using $(nproc) cores)..."
make -j$(nproc)

cd - > /dev/null

# Verify
BINDING_FILE=$(find crm_ml_rl -name "crm_cpp*.so" 2>/dev/null | head -1)

if [ -z "$BINDING_FILE" ]; then
    print_error "Build failed - binding not found"
    exit 1
fi

print_status "C++ bindings built: $BINDING_FILE"

# Test import
print_info "Testing import..."
if python3 -c "import crm_cpp; print('Module version:', crm_cpp.__version__)" 2>/dev/null; then
    print_status "Import successful!"
else
    print_error "Import failed - check dependencies"
    exit 1
fi

echo ""
echo "=========================================="
echo -e "${GREEN}Build Complete!${NC}"
echo "=========================================="
echo ""
echo "Usage:"
echo "  python3 -c 'import crm_cpp'"
echo "  from crm_ml_rl import ForwardKinematicsAPI"
echo ""
echo "Options:"
echo "  ./scripts/build_cpp.sh          # Incremental build"
echo "  ./scripts/build_cpp.sh --clean  # Clean + rebuild"
