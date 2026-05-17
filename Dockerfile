FROM debian:bookworm-slim AS builder

RUN apt-get update && apt-get install -y \
    g++ curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN c++ -std=c++20 -O2 -I. \
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
    server/server.cpp \
    -o cfd_server -lpthread

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y curl ca-certificates && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY --from=builder /app/cfd_server .
COPY --from=builder /app/webui ./webui

ENV PORT=8000
EXPOSE 8000

CMD ["./cfd_server"]
