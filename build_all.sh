#!/bin/bash
set -e

cd /workspaces/CRM_Dynamics

echo "=========================================="
echo "FULL BUILD PROCESS"
echo "=========================================="

# Clean
echo ""
echo "[1/7] Cleaning previous builds..."
rm -rf build_debug build_release
mkdir -p build_debug build_release

# Configure Debug
echo ""
echo "[2/7] Configuring CMake (Debug)..."
cd build_debug
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_STANDARD=17 \
      -G Ninja \
      -DPython3_EXECUTABLE=$(which python3) \
      .. || exit 1
cd ..

# Configure Release
echo ""
echo "[3/7] Configuring CMake (Release)..."
cd build_release
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=17 \
      -G Ninja \
      -DPython3_EXECUTABLE=$(which python3) \
      .. || exit 1
cd ..

# Build Debug
echo ""
echo "[4/7] Building Debug..."
cd build_debug
ninja
cd ..

# Build Python module
echo ""
echo "[5/7] Building Python Module..."
cd build_debug
ninja crm_cpp || echo "Warning: crm_cpp build skipped"
cd ..

# Test executables
echo ""
echo "[6/7] Testing Executables..."
./build_debug/CRMTest > /dev/null && echo "  ✓ CRMTest works" || echo "  ✗ CRMTest failed"
./build_debug/CRMDYNTest > /dev/null && echo "  ✓ CRMDYNTest works" || echo "  ✗ CRMDYNTest failed"

# Test Python
echo ""
echo "[7/7] Testing Python Module..."
source .venv/bin/activate
python3 << 'EOF'
import sys
sys.path.insert(0, './build_debug')
try:
    import crm_cpp
    print("  ✓ crm_cpp module works")
except:
    print("  ✗ crm_cpp module not found")
EOF

echo ""
echo "=========================================="
echo "✓ BUILD COMPLETE!"
echo "=========================================="
echo ""
echo "Quick test commands:"
echo "  ./build_debug/CRMTest"
echo "  ./build_debug/CRMDYNTest"
echo "  python3 -c 'import sys; sys.path.insert(0, \"./build_debug\"); import crm_cpp; print(crm_cpp.version_info())'"
