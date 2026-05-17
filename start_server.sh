#!/usr/bin/env bash
set -e
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo ""
echo "  Discreetize CFD Mesh — Web / Production Build"
echo "  =============================================="
echo "  Deploy: Railway (Dockerfile) + Lovable + Supabase"
echo ""

echo "==> Compiling..."
c++ -std=c++20 -O2 -I. \
    geometry/geometry_engine.cpp \
    geometry/bvh.cpp \
    geometry/primitives.cpp \
    geometry/mesh_sdf.cpp \
    geometry/bifilar.cpp \
    meshing/mesh_engine.cpp \
    meshing/quad_mesh.cpp \
    analysis/curvature.cpp \
    slicer/slicer.cpp \
    export/mesh_exporter.cpp \
    supabase/supabase_client.cpp \
    cfd/cfd_engine.cpp \
    shaders/shader_library.cpp \
    shaders/node_compiler.cpp \
    render/software_renderer.cpp \
    server/server.cpp \
    -o cfd_server -lpthread

# Railway injects PORT env var; argv[1] used for local testing
PORT="${PORT:-${1:-8000}}"
echo "==> Starting server on port $PORT"
./cfd_server "$PORT"
