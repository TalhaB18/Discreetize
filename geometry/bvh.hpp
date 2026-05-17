#pragma once
#include "types.hpp"
#include <vector>
#include <cstdint>

namespace cfd {

struct AABB {
    Vec3d lo, hi;
    bool intersectsZSlab(double z_lo, double z_hi) const;
    void expand(const Vec3d& p);
    static AABB merge(const AABB& a, const AABB& b);
    static AABB fromTriangle(const Vec3d& a, const Vec3d& b, const Vec3d& c);
};

class BVH {
public:
    void build(const SurfaceManifold& m);

    // Returns indices of triangles whose AABB overlaps the Z slab [z-eps, z+eps].
    // eps defaults to 0 — exact slab query.
    std::vector<uint32_t> queryZ(double z, double eps = 0.0) const;

    // Returns the index of the nearest triangle to point p, and sets dist_out
    // to the unsigned distance from p to that triangle's surface.
    uint32_t nearestTriangle(const Vec3d& p, double& dist_out) const;

    bool empty() const { return nodes_.empty(); }

private:
    struct Node {
        AABB  bounds;
        int   left  = -1;   // -1 = leaf
        int   right = -1;
        uint32_t face_begin = 0;  // leaf: index into leaf_faces_
        uint32_t face_count = 0;
    };

    std::vector<Node>     nodes_;
    std::vector<uint32_t> leaf_faces_;  // face indices stored contiguously for all leaves
    const SurfaceManifold* mesh_ = nullptr;

    // Build a node covering face_ids[begin..end). Returns node index.
    int buildNode(std::vector<uint32_t>& face_ids, int begin, int end);

    static constexpr int LEAF_SIZE = 8;  // max faces per leaf
};

} // namespace cfd
