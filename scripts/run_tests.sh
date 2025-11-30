#!/bin/bash
set -e

echo "=========================================="
echo "Running CRM Catheter Test Suite"
echo "=========================================="

echo ""
echo "Running C++ tests..."
cd build
ctest --output-on-failure --verbose
cd ..

echo ""
echo "Running Python tests..."
source .venv/bin/activate
pytest tests/python/ -v --cov=crm_catheter --cov-report=html --cov-report=term

echo ""
echo "=========================================="
echo "✅ All tests passed!"
echo "=========================================="
