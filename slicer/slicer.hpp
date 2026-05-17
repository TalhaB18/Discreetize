#pragma once
#include "../geometry/types.hpp"
#include <vector>
#include <array>
#include <string>

namespace cfd {

using Vec2d = std::array<double, 2>;

struct Contour {
    std::vector<Vec2d> points;   // ordered closed loop in XY plane
    double z;
};

struct SliceResult {
    std::vector<Contour> contours;
    double z_level;
};

class Slicer {
public:
    SliceResult sliceAtZ(const SurfaceManifold& m, double z) const;

    // z_min/z_max auto-detected from mesh AABB when both are 0.0
    std::vector<SliceResult> sliceRange(const SurfaceManifold& m,
                                         double z_min, double z_max,
                                         double z_interval) const;

    void exportSVG(const std::vector<SliceResult>& slices,
                   const std::string& path) const;

    void exportPolylines(const std::vector<SliceResult>& slices,
                          const std::string& path) const;
};

} // namespace cfd
