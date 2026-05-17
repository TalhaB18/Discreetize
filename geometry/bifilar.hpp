#pragma once
#include "types.hpp"
#include <functional>

namespace cfd {

struct BifilarConfig {
    double R           = 1.0;   // outer cylinder radius
    double r           = 0.44;  // inner cylinder radius
    double offset      = 0.48;  // inner circle centre offset from axis
    double H           = 4.0;   // total height
    double revolutions = 2.0;   // full rotations over H
    int    resolution  = 50;    // MC voxels per unit length
};

// SDF for the fluid domain (negative = inside fluid):
//   max(sdf_outer, -sdf_inner1, -sdf_inner2, -z, z-H)
std::function<double(double,double,double)> makeBifilarSDF(const BifilarConfig& cfg);

// Generate watertight surface mesh of the full fluid domain via Marching Cubes.
SurfaceManifold generateBifilarSurface(const BifilarConfig& cfg);

// Generate one solid inner helical cylinder (which = 1 or 2) for visualisation.
// The two cylinders are 180° apart and wind around the axis as height increases.
SurfaceManifold generateBifilarInnerCylinder(const BifilarConfig& cfg, int which);

} // namespace cfd
