#!/bin/bash
# Quick build script for AstDyn C++

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}╔════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║    AstDyn C++ - Quick Build Script        ║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════════╝${NC}"
echo ""

# Parse arguments
BUILD_TYPE="Release"
BUILD_TESTS="ON"
BUILD_EXAMPLES="ON"
BUILD_DOCS="OFF"
CLEAN_BUILD=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--debug)
            BUILD_TYPE="Debug"
            shift
            ;;
        -r|--release)
            BUILD_TYPE="Release"
            shift
            ;;
        --no-tests)
            BUILD_TESTS="OFF"
            shift
            ;;
        --no-examples)
            BUILD_EXAMPLES="OFF"
            shift
            ;;
        --docs)
            BUILD_DOCS="ON"
            shift
            ;;
        --clean)
            CLEAN_BUILD=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -d, --debug         Build in Debug mode (default: Release)"
            echo "  -r, --release       Build in Release mode"
            echo "  --no-tests          Don't build tests"
            echo "  --no-examples       Don't build examples"
            echo "  --docs              Build documentation"
            echo "  --clean             Clean build directory first"
            echo "  -h, --help          Show this help message"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

# Clean build directory if requested
if [ "$CLEAN_BUILD" = true ]; then
    echo -e "${YELLOW}Cleaning build directory...${NC}"
    rm -rf build
fi

# Create build directory
mkdir -p build
cd build

# Configure
echo -e "${YELLOW}Configuring CMake...${NC}"
cmake .. \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DASTDYN_BUILD_TESTS=$BUILD_TESTS \
    -DASTDYN_BUILD_EXAMPLES=$BUILD_EXAMPLES \
    -DASTDYN_BUILD_DOCS=$BUILD_DOCS \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

# Build
echo ""
echo -e "${YELLOW}Building...${NC}"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
cmake --build . -j$NPROC

# Run tests if enabled
if [ "$BUILD_TESTS" = "ON" ]; then
    echo ""
    echo -e "${YELLOW}Running tests...${NC}"
    ctest --output-on-failure
fi

# Build docs if enabled
if [ "$BUILD_DOCS" = "ON" ]; then
    echo ""
    echo -e "${YELLOW}Building documentation...${NC}"
    make docs
fi

echo ""
echo -e "${GREEN}╔════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║         Build completed successfully!      ║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════════╝${NC}"
echo ""
echo -e "Build type:    ${GREEN}$BUILD_TYPE${NC}"
echo -e "Tests built:   ${GREEN}$BUILD_TESTS${NC}"
echo -e "Examples built: ${GREEN}$BUILD_EXAMPLES${NC}"
echo -e "Docs built:    ${GREEN}$BUILD_DOCS${NC}"
echo ""

if [ "$BUILD_EXAMPLES" = "ON" ]; then
    echo -e "${YELLOW}Run example:${NC} ./examples/example_basic"
fi

if [ "$BUILD_TESTS" = "ON" ]; then
    echo -e "${YELLOW}Run tests:${NC} ctest --output-on-failure"
fi

echo ""
