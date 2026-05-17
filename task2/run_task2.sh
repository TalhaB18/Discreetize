#!/bin/bash
# Build and run Task 2 — Twisted Bifilar Solid pipeline
set -e
cd "$(dirname "$0")/.."

BUILD_DIR="build"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON 2>&1 | tail -5
cmake --build "$BUILD_DIR" --target task2_pipeline -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"

echo ""
echo "Running task2_pipeline..."
"$BUILD_DIR/task2_pipeline" --output-dir output/task2 "$@"

echo ""
echo "Generating cross-section visualisation..."
python3 task2/visualize_paraview.py --matplotlib-only

echo ""
echo "Done. Open output/task2/ to view results."
echo "For full ParaView renders: pvpython task2/visualize_paraview.py"
