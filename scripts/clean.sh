#!/bin/bash
# Clean up generated files

echo "Cleaning up..."

# Remove data
rm -rf data/*.npz

# Remove models
rm -rf models/*.pth

# Remove logs
rm -rf logs/*

# Remove results
rm -rf results/*

# Remove build artifacts
rm -rf crm_ml_rl/bindings/build/
find crm_ml_rl -name "crm_cpp*.so" -delete

# Remove Python cache
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -type f -name "*.pyc" -delete

# Remove test artifacts
rm -rf .pytest_cache/
rm -rf htmlcov/
rm -f .coverage

# Remove generated images
rm -f *.png

echo "✓ Cleanup complete!"
