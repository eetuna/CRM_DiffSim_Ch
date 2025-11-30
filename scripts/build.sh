#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

BUILD_TYPE="Release"
PYTHON_BINDINGS=ON
TESTING=ON
EXAMPLES=ON
CLEAN=false
INSTALL=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --debug) BUILD_TYPE="Debug"; shift ;;
        --release) BUILD_TYPE="Release"; shift ;;
        --no-python) PYTHON_BINDINGS=OFF; shift ;;
        --no-tests) TESTING=OFF; shift ;;
        --no-examples) EXAMPLES=OFF; shift ;;
        --clean) CLEAN=true; shift ;;
        --install) INSTALL=true; shift ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  --debug          Build in Debug mode"
            echo "  --release        Build in Release mode (default)"
            echo "  --no-python      Disable Python bindings"
            echo "  --no-tests       Disable tests"
            echo "  --no-examples    Disable examples"
            echo "  --clean          Clean before building"
            echo "  --install        Install after building"
            exit 0
            ;;
        *) echo -e "${RED}Unknown option: $1${NC}"; exit 1 ;;
    esac
done

cd "$PROJECT_ROOT"

if [ "$CLEAN" = true ]; then
    echo -e "${YELLOW}Cleaning build directories...${NC}"
    rm -rf build build_debug build_release
fi

BUILD_DIR="build_$(echo $BUILD_TYPE | tr '[:upper:]' '[:lower:]')"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

echo -e "${GREEN}Configuring CMake ($BUILD_TYPE)...${NC}"
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
      -DBUILD_PYTHON_BINDINGS=$PYTHON_BINDINGS \
      -DBUILD_TESTING=$TESTING \
      -DBUILD_EXAMPLES=$EXAMPLES \
      -G Ninja \
      ..

echo -e "${GREEN}Building...${NC}"
ninja

if [ "$TESTING" = "ON" ]; then
    echo -e "${GREEN}Running tests...${NC}"
    ctest --output-on-failure
fi

if [ "$INSTALL" = true ]; then
    echo -e "${GREEN}Installing...${NC}"
    sudo ninja install
fi

echo -e "${GREEN}Build complete!${NC}"
