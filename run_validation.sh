#!/bin/bash
# run_validation.sh

echo "=========================================="
echo "  ITALOccultLibrary Validation Pipeline"
echo "=========================================="

# 1. Fetch Fresh Data from JPL Horizons
echo "[1/3] Fetching Ground Truth Data from JPL..."
python3 fetch_jpl_horizons.py
if [ $? -ne 0 ]; then
    echo "❌ Failed to fetch data."
    exit 1
fi

# 2. Compile Test
echo "[2/3] Compiling Validation Test..."
g++ -std=c++17 -O2 -I./italoccultlibrary/include -I./include -I./astdyn/include \
    -I/opt/homebrew/include -I/usr/local/include/eigen3 -I/opt/homebrew/include/eigen3 \
    test_jpl_validation.cpp italoccultlibrary/src/astdyn_wrapper.cpp \
    $(find astdyn/src -name "*.cpp" | grep -v "main.cpp" | grep -v "astdyn_convert.cpp") \
    -o test_jpl_validation -lstdc++

if [ $? -ne 0 ]; then
    echo "❌ Compilation failed."
    exit 1
fi

# 3. Run Test
echo "[3/3] Running Validation Test..."
./test_jpl_validation

RET=$?
if [ $RET -eq 0 ]; then
    echo ""
    echo "✅ PIPELINE SUCCESS: Library validated against JPL Horizons."
else
    echo ""
    echo "❌ PIPELINE FAILED: Accuracy requirements not met."
fi

exit $RET
