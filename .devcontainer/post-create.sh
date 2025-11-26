#!/bin/bash
set -e

echo "=========================================="
echo "Post-Create Script: Setting up environment"
echo "WSL 2 / Ubuntu 22.04 Edition"
echo "=========================================="

# ============================================================
# FIX FOR WSL 2: Update package list (may be stale)
# ============================================================
echo "Updating package list (WSL 2)..."
sudo apt-get update --allow-insecure-repositories || true

# ============================================================
# Create Python virtual environment in workspace
# ============================================================
if [ ! -d ".venv" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv .venv
fi

# Activate venv and install dependencies
source .venv/bin/activate

echo "Installing Python dependencies from requirements.txt (if exists)..."
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    pip install --upgrade pip setuptools wheel
    pip install \
        numpy \
        scipy \
        matplotlib \
        scikit-learn \
        pandas \
        torch \
        pybind11[global]
fi

# ============================================================
# Create build directory structure
# ============================================================
echo "Setting up C++ build directories..."
mkdir -p build_debug build_release

# ============================================================
# CMake Configuration for Debug Build
# ============================================================
echo "Configuring CMake (Debug)..."
cd build_debug
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_STANDARD=17 \
      -G Ninja \
      ..
cd ..

# ============================================================
# CMake Configuration for Release Build
# ============================================================
echo "Configuring CMake (Release)..."
cd build_release
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -G Ninja \
      ..
cd ..

# ============================================================
# WSL 2 SPECIFIC: Fix line endings (CRLF → LF)
# ============================================================
echo "Normalizing line endings for WSL 2..."
find . -name "*.cpp" -o -name "*.hpp" -o -name "*.h" | xargs dos2unix 2>/dev/null || true

# ============================================================
# EXPLICITLY BUILD PYTHON BINDINGS
# ============================================================
echo ""
echo "Building Python bindings (crm_cpp)..."
cd build_debug
ninja crm_cpp 2>/dev/null || echo "Warning: crm_cpp build skipped (pybind11 or source file may be missing)"
cd ..

echo "=========================================="
echo "✓ Development environment ready!"
echo "=========================================="
echo ""
echo "Quick commands:"
echo "  Build (Debug):    cd build_debug && ninja"
echo "  Build (Release):  cd build_release && ninja"
echo "  Run tests:        ./build_debug/CRMTest"
echo "  Run dynamics:     ./build_debug/CRMDYNTest"
echo ""
echo "Python (in .venv):"
echo "  source .venv/bin/activate"
echo "  python your_script.py"
echo ""
echo "Test Python bindings:"
echo "  source .venv/bin/activate"
echo "  python -c 'import sys; sys.path.insert(0, \"./build_debug\"); import crm_cpp; print(crm_cpp.version_info())'"
echo ""
echo "To enter container:"
echo "  Ctrl+Shift+P → Dev Containers: Open Folder in Container"
echo ""