#pragma once
#include "types.hpp"
#include "bvh.hpp"
#include <functional>

namespace cfd {

// Signed Distance Field evaluated from an arbitrary closed SurfaceManifold.
// Positive outside, negative inside.
class MeshSDF {
public:
    explicit MeshSDF(const SurfaceManifold& mesh);

    double eval(const Vec3d& p) const;

    std::function<double(double,double,double)> asLambda() const {
        return [this](double x, double y, double z){ return eval({x,y,z}); };
    }

    void bounds(Vec3d& lo, Vec3d& hi) const { lo = bbox_lo_; hi = bbox_hi_; }

private:
    const SurfaceManifold& mesh_;
    BVH bvh_;
    Vec3d bbox_lo_, bbox_hi_;

    bool isInside(const Vec3d& p) const;
    int  raycastCount(const Vec3d& p, const Vec3d& dir) const;
};

} // namespace cfd
