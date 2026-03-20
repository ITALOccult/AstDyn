#!/bin/bash
# Quick compile test for AstDyn manual
# Usage: ./compile_test.sh [en|it]

LANG=${1:-en}

if [ "$LANG" = "en" ]; then
    echo "Testing English manual compilation..."
    cd en/
elif [ "$LANG" = "it" ]; then
    echo "Testing Italian manual compilation..."
    cd it/
else
    echo "Usage: $0 [en|it]"
    exit 1
fi

# Quick single-pass compilation
make quick

if [ $? -eq 0 ]; then
    echo "✅ Compilation successful!"
    echo "Output: ../../build/manual/main*.pdf"
    ls -lh ../../build/manual/*.pdf 2>/dev/null
else
    echo "❌ Compilation failed. Check errors above."
    exit 1
fi
