#!/bin/bash
# Run all tests

set -e

echo "=========================================="
echo "Running Tests"
echo "=========================================="

## Run pytest with coverage
# pytest tests/ -v --cov=crm_ml_rl --cov-report=term --cov-report=html
# Run pytest with coverage (use `python -m pytest` for interpreter-resolved invocation)
python -m pytest tests/ -v --cov=crm_ml_rl --cov-report=term --cov-report=html


echo ""
echo "✓ All tests passed!"
echo "  Coverage report: htmlcov/index.html"
