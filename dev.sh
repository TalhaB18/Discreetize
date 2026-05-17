#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# dev.sh — Hot-reload development server
#
# • Compiles and starts cfd_server
# • Watches *.cpp / *.hpp for changes → auto-recompiles and restarts
# • The C++ server itself watches webui/ → sends SSE to browser → auto-reload
#
# Usage:  ./dev.sh [port]   (default port 8000)
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

PORT="${1:-8000}"
SERVER_PID=""

BUILD_CMD=(
    c++ -std=c++20 -O2 -I.
    geometry/geometry_engine.cpp
    geometry/bvh.cpp
    geometry/primitives.cpp
    geometry/mesh_sdf.cpp
    geometry/bifilar.cpp
    meshing/mesh_engine.cpp
    meshing/quad_mesh.cpp
    analysis/curvature.cpp
    slicer/slicer.cpp
    export/mesh_exporter.cpp
    supabase/supabase_client.cpp
    cfd/cfd_engine.cpp
    cfd/lbm_sim.cpp
    cfd/flow3d.cpp
    shaders/shader_library.cpp
    shaders/node_compiler.cpp
    render/software_renderer.cpp
    server/server.cpp
    -o cfd_server -lpthread -lcurl
)

build() {
    echo "==> Building..."
    if "${BUILD_CMD[@]}" 2>&1; then
        echo "==> Build OK"
        return 0
    else
        echo "==> Build FAILED — server not restarted"
        return 1
    fi
}

start_server() {
    ./cfd_server "$PORT" &
    SERVER_PID=$!
    echo "==> Server started  pid=$SERVER_PID  http://localhost:$PORT"
}

stop_server() {
    if [ -n "$SERVER_PID" ] && kill -0 "$SERVER_PID" 2>/dev/null; then
        kill "$SERVER_PID" 2>/dev/null
        wait "$SERVER_PID" 2>/dev/null || true
        SERVER_PID=""
    fi
}

# Snapshot mtimes of all watched C++/H files into a temp file
STAMP_FILE="$(mktemp /tmp/cfd_dev_stamp.XXXXXX)"
snapshot() {
    find . \( -name "*.cpp" -o -name "*.hpp" \) \
        ! -path "./webui/*" ! -path "./.git/*" \
        -exec stat -f "%m %N" {} \; 2>/dev/null | sort > "$STAMP_FILE"
}

changed() {
    find . \( -name "*.cpp" -o -name "*.hpp" \) \
        ! -path "./webui/*" ! -path "./.git/*" \
        -exec stat -f "%m %N" {} \; 2>/dev/null | sort | diff - "$STAMP_FILE" > /dev/null 2>&1
    # diff returns 0 if same, 1 if different — invert for "changed"
    [ $? -ne 0 ]
}

cleanup() {
    stop_server
    rm -f "$STAMP_FILE"
    echo ""
    echo "==> Dev server stopped"
    exit 0
}
trap cleanup INT TERM

echo "==> CFD Mesh Dev Server"
echo "==> Watching *.cpp / *.hpp for changes (webui/ watched by C++ server)"

build && start_server
snapshot

while true; do
    sleep 1
    if changed; then
        echo "==> Source changed — rebuilding..."
        stop_server
        if build; then
            start_server
        fi
        snapshot
    fi
done
