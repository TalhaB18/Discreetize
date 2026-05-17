#pragma once
#include "../geometry/types.hpp"
#include <vector>
#include <array>

namespace cfd {

struct CurvatureResult {
    std::vector<double>                 mean_curvature; // per vertex, signed H
    std::vector<std::array<uint8_t,3>> jet_colors;     // RGB [0,255] per vertex
};

class CurvatureAnalyzer {
public:
    CurvatureResult compute(const SurfaceManifold& m) const;
private:
    static double cotWeight(const Vec3d& p, const Vec3d& q, const Vec3d& r);
    static std::array<uint8_t,3> jetColor(double t);
};

} // namespace cfd
