cmake_minimum_required(VERSION 3.20)
project(CFDMeshWrapper LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-O2 -Wall -Wextra -march=native)
endif()

# Library target
add_library(cfd_mesh STATIC
    geometry/bvh.cpp
    geometry/geometry_engine.cpp
    geometry/primitives.cpp
    geometry/mesh_sdf.cpp
    geometry/bifilar.cpp
    meshing/mesh_engine.cpp
    meshing/quad_mesh.cpp
    cfd/cfd_engine.cpp
    analysis/curvature.cpp
    slicer/slicer.cpp
    export/mesh_exporter.cpp
    supabase/supabase_client.cpp
    shaders/shader_library.cpp
    render/software_renderer.cpp
)
target_include_directories(cfd_mesh PUBLIC ${CMAKE_SOURCE_DIR})

# Main CLI executable
add_executable(cfd_pipeline main.cpp)
target_link_libraries(cfd_pipeline PRIVATE cfd_mesh)

# Task 1 — Boolean Subtraction + Tetrahedral Meshing
add_executable(task1_pipeline task1/task1_pipeline.cpp)
target_link_libraries(task1_pipeline PRIVATE cfd_mesh)

# Task 2 — Twisted Bifilar Solid
add_executable(task2_pipeline task2/task2_pipeline.cpp)
target_link_libraries(task2_pipeline PRIVATE cfd_mesh)

# HTTP server
add_executable(cfd_server server/server.cpp)
target_link_libraries(cfd_server PRIVATE cfd_mesh)
target_compile_options(cfd_server PRIVATE -pthread)
target_link_options(cfd_server PRIVATE -pthread)

# Optional: find Gmsh
find_library(GMSH_LIB gmsh HINTS /usr/local/lib /opt/homebrew/lib)
if(GMSH_LIB)
    message(STATUS "Gmsh found: ${GMSH_LIB}")
    target_link_libraries(cfd_mesh PUBLIC ${GMSH_LIB})
    target_compile_definitions(cfd_mesh PUBLIC HAS_GMSH=1)
endif()
