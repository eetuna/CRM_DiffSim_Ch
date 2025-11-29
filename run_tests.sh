#!/bin/bash
set -e

echo "=========================================="
echo "COMPREHENSIVE TEST SUITE"
echo "=========================================="

# Test 1: C++ Build
echo ""
echo "[1/5] Testing C++ Kinematics..."
cd build_debug && ninja && cd ..
./build_debug/CRMTest

# Test 2: C++ Dynamics
echo ""
echo "[2/5] Testing C++ Dynamics..."
cd build_debug && ninja CRMDYNTest && cd ..
./build_debug/CRMDYNTest

# Test 3: Python Bindings
echo ""
echo "[3/5] Testing Python Bindings..."
source .venv/bin/activate
python3 << 'EOF'
import sys
sys.path.insert(0, './build_debug')
import crm_cpp
print("✓ Python bindings work!")
EOF

# Test 4: ML Models
echo ""
echo "[4/5] Testing ML Models..."
python3 test_ml_models.py

# Test 5: Integration
echo ""
echo "[5/5] Testing Full Integration..."
python3 test_integration.py

echo ""
echo "=========================================="
echo "✓ ALL TESTS PASSED!"
echo "=========================================="
