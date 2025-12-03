#!/bin/bash
# Runs automatically after DevContainer is created

set -e  # Exit on error
set -x  # Show each command (DEBUG MODE)


echo "=========================================="
echo "DevContainer Post-Create Setup"
echo "=========================================="
echo ""

# Install Python packages
echo "→ Installing Python packages..."
REQUIREMENTS_FILE="/workspaces/CRM_Dynamics/requirements.txt"
if [ -f "$REQUIREMENTS_FILE" ]; then
    #pip install --upgrade pip
    pip install -r requirements.txt
    echo "✓ Python packages installed"
else
    echo "⚠ requirements.txt not found"
fi
echo ""

# Make scripts executable
echo "→ Making scripts executable..."
if [ -d "scripts" ]; then
    chmod +x scripts/*.sh
    echo "✓ Scripts are executable"
else
    echo "⚠ scripts/ directory not found"
fi
echo ""

# Verify system tools
echo "→ Verifying installed tools..."
echo "  Python: $(python --version 2>&1)"
echo "  CMake: $(cmake --version | head -1)"
echo "  g++: $(g++ --version | head -1)"
echo "  Eigen3: $(pkg-config --modversion eigen3 2>/dev/null || echo 'installed')"
echo ""

# Check CRM_PATH
echo "→ Checking CRM configuration..."
if [ -z "$CRM_PATH" ]; then
    echo "⚠ CRM_PATH not set"
else
    echo "✓ CRM_PATH: $CRM_PATH"
    if [ -d "$CRM_PATH" ]; then
        echo "✓ CRM directory exists"
    else
        echo "⚠ CRM directory not found at: $CRM_PATH"
    fi
fi
echo ""

echo "=========================================="
echo "✓ Post-Create Setup Complete!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  Run: make help"
echo "  Run: make setup"
echo ""
