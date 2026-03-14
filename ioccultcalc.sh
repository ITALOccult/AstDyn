#!/bin/bash
# ioccultcalc wrapper script

# Get the directory where the script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
BUILD_DIR="${SCRIPT_DIR}/build"

# Build if binary doesn't exist
if [ ! -f "${BUILD_DIR}/astdyn/tools/bin/ioccultcalc" ]; then
    echo "Building ioccultcalc..."
    cmake -B "${BUILD_DIR}" "${SCRIPT_DIR}" && cmake --build "${BUILD_DIR}" --target ioccultcalc -j$(sysctl -n hw.ncpu)
fi

# Run the tool
"${BUILD_DIR}/astdyn/tools/bin/ioccultcalc" "$@"
