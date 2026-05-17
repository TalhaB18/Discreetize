cmake_minimum_required(VERSION 3.20)
project(cfd_geometry CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Geometry library
add_library(cfd_geometry STATIC
    geometry_engine.cpp
)
target_include_directories(cfd_geometry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Optional: enable stricter warnings
target_compile_options(cfd_geometry PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
        -Wall -Wextra -Wpedantic -Wshadow>
)

# Unit test executable
add_executable(geometry_test geometry_test.cpp)
target_link_libraries(geometry_test PRIVATE cfd_geometry)
