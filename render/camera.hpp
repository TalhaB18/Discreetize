#pragma once
#include <cmath>

namespace render {

struct Camera {
    double eye[3]    = {2.0, 2.0, 2.0};
    double target[3] = {0.0, 0.0, 0.0};
    double up[3]     = {0.0, 1.0, 0.0};
    double fov_deg   = 45.0;
    int    width     = 800;
    int    height    = 600;
};

} // namespace render
