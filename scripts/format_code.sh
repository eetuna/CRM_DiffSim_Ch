#!/bin/bash
set -e

echo "Formatting C++ code with clang-format..."
find src numerical main python -name "*.cpp" -o -name "*.hpp" -o -name "*.h" | \
    xargs clang-format -i

echo "Formatting Python code with black..."
source .venv/bin/activate
black python/ tests/ examples/

echo "Sorting Python imports with isort..."
isort python/ tests/ examples/

echo "✅ Code formatting complete!"
