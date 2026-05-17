#pragma once
#include <vector>
#include <cstdint>
#include "camera.hpp"

namespace render {

// Render a triangle soup (N*9 floats: 3 verts × 3 coords per tri) with Phong
// shading into a PNG image returned as raw bytes.
//
// bg_color / mesh_color are 0xRRGGBB packed uint32_t values.
std::vector<uint8_t> render_stl(
    const std::vector<float>& triangles,
    const Camera& cam,
    uint32_t bg_color   = 0x1a1a2e,
    uint32_t mesh_color = 0x6b58cc
);

} // namespace render
