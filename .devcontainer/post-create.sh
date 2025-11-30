#!/bin/bash
set -e

echo "========================================"
echo "Post-Create: CRM Catheter Development Setup"
echo "========================================"

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Ensure we're in the workspace directory
cd "${WORKSPACE_FOLDER:-/workspace}" || cd /workspace || { echo "ERROR: Cannot find workspace"; exit 1; }

echo ""
echo "Step 1: Python Virtual Environment"
echo "------------------------------------"

# Check if Python is available
if ! command_exists python3; then
    echo "ERROR: Python3 not found. Install it first."
    exit 1
fi

# Create venv if it doesn't exist
if [ ! -d ".venv" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv .venv
else
    echo "Virtual environment already exists"
fi

# Activate venv
if [ -f ".venv/bin/activate" ]; then
    source .venv/bin/activate
else
    echo "ERROR: Cannot find .venv/bin/activate"
    exit 1
fi

echo ""
echo "Step 2: Install Python Dependencies"
echo "------------------------------------"

# Upgrade pip
pip install --upgrade pip setuptools wheel

# Install dependencies if requirements files exist
if [ -f "requirements.txt" ]; then
    echo "Installing requirements.txt..."
    pip install -r requirements.txt
else
    echo "WARNING: requirements.txt not found, skipping"
fi

if [ -f "requirements-dev.txt" ]; then
    echo "Installing requirements-dev.txt..."
    pip install -r requirements-dev.txt
else
    echo "WARNING: requirements-dev.txt not found, skipping"
fi

echo ""
echo "Step 3: CMake Configuration"
echo "------------------------------------"

# Check if CMake is available
if ! command_exists cmake; then
    echo "ERROR: CMake not found. Install it first."
    exit 1
fi

# Create build directory
mkdir -p build
cd build

# Configure with CMake
echo "Configuring CMake..."
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_PYTHON_BINDINGS=ON \
    -DBUILD_TESTING=ON \
    -DBUILD_EXAMPLES=ON \
    -DENABLE_OPENMP=ON \
    -G Ninja \
    .. || {
    echo "WARNING: CMake configuration failed, continuing anyway..."
    cd ..
}

cd ..

echo ""
echo "Step 4: Build C++ Library"
echo "------------------------------------"

if [ -d "build" ] && command_exists ninja; then
    cd build
    echo "Building with Ninja..."
    ninja || {
        echo "WARNING: Build failed, but continuing..."
    }
    cd ..
else
    echo "WARNING: Skipping build (ninja or build dir not found)"
fi

echo ""
echo "Step 5: Install Python Package"
echo "------------------------------------"

# Install package in development mode
if [ -f "setup.py" ]; then
    echo "Installing Python package in development mode..."
    pip install -e . || {
        echo "WARNING: Python package installation failed"
    }
else
    echo "WARNING: setup.py not found, skipping package install"
fi

echo ""
echo "Step 6: Git Configuration"
echo "------------------------------------"

# Configure git safe directory
if command_exists git; then
    git config --global --add safe.directory "${PWD}"
    echo "Git safe directory configured"
else
    echo "WARNING: Git not found"
fi

echo ""
echo "========================================"
echo "✓ Setup Complete!"
echo "========================================"
echo ""
echo "Quick Commands:"
echo "  Build C++:      cd build && ninja"
echo "  Run tests:      cd build && ctest"
echo "  Python import:  python -c 'import crm_cpp'"
echo "  Activate venv:  source .venv/bin/activate"
echo ""
