#include "bvh.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>

namespace cfd {

// ---------------------------------------------------------------------------
// AABB
// ---------------------------------------------------------------------------

bool AABB::intersectsZSlab(double z_lo, double z_hi) const {
    return lo[2] <= z_hi && hi[2] >= z_lo;
}

void AABB::expand(const Vec3d& p) {
    for (int k = 0; k < 3; ++k) {
        if (p[k] < lo[k]) lo[k] = p[k];
        if (p[k] > hi[k]) hi[k] = p[k];
    }
}

AABB AABB::merge(const AABB& a, const AABB& b) {
    AABB result;
    for (int k = 0; k < 3; ++k) {
        result.lo[k] = std::min(a.lo[k], b.lo[k]);
        result.hi[k] = std::max(a.hi[k], b.hi[k]);
    }
    return result;
}

AABB AABB::fromTriangle(const Vec3d& a, const Vec3d& b, const Vec3d& c) {
    AABB box;
    for (int k = 0; k < 3; ++k) {
        box.lo[k] = std::min({a[k], b[k], c[k]});
        box.hi[k] = std::max({a[k], b[k], c[k]});
    }
    return box;
}

// ---------------------------------------------------------------------------
// BVH
// ---------------------------------------------------------------------------

void BVH::build(const SurfaceManifold& m) {
    mesh_ = &m;
    nodes_.clear();
    leaf_faces_.clear();

    const uint32_t nf = static_cast<uint32_t>(m.faces.size());
    if (nf == 0) return;

    std::vector<uint32_t> face_ids(nf);
    for (uint32_t i = 0; i < nf; ++i) face_ids[i] = i;

    buildNode(face_ids, 0, static_cast<int>(nf));
}

int BVH::buildNode(std::vector<uint32_t>& face_ids, int begin, int end) {
    // Compute AABB for all faces in [begin, end).
    const double kInf = std::numeric_limits<double>::infinity();
    AABB bounds;
    bounds.lo = { kInf,  kInf,  kInf};
    bounds.hi = {-kInf, -kInf, -kInf};

    for (int i = begin; i < end; ++i) {
        const auto& tri = mesh_->faces[face_ids[i]];
        const Vec3d& p0 = mesh_->vertices[tri.vi[0]].pos;
        const Vec3d& p1 = mesh_->vertices[tri.vi[1]].pos;
        const Vec3d& p2 = mesh_->vertices[tri.vi[2]].pos;
        bounds.expand(p0);
        bounds.expand(p1);
        bounds.expand(p2);
    }

    // Allocate node slot now (before recursion may reallocate nodes_).
    const int node_idx = static_cast<int>(nodes_.size());
    nodes_.emplace_back();
    nodes_[node_idx].bounds = bounds;

    if (end - begin <= LEAF_SIZE) {
        // Leaf node.
        nodes_[node_idx].left       = -1;
        nodes_[node_idx].right      = -1;
        nodes_[node_idx].face_begin = static_cast<uint32_t>(leaf_faces_.size());
        nodes_[node_idx].face_count = static_cast<uint32_t>(end - begin);

        for (int i = begin; i < end; ++i) {
            leaf_faces_.push_back(face_ids[i]);
        }
        return node_idx;
    }

    // Find longest axis.
    Vec3d extent = {
        bounds.hi[0] - bounds.lo[0],
        bounds.hi[1] - bounds.lo[1],
        bounds.hi[2] - bounds.lo[2]
    };
    int axis = 0;
    if (extent[1] > extent[axis]) axis = 1;
    if (extent[2] > extent[axis]) axis = 2;

    // Sort faces by centroid along chosen axis.
    std::sort(face_ids.begin() + begin, face_ids.begin() + end,
        [&](uint32_t fa, uint32_t fb) {
            const auto& ta = mesh_->faces[fa];
            const auto& tb = mesh_->faces[fb];
            double ca = (mesh_->vertices[ta.vi[0]].pos[axis] +
                         mesh_->vertices[ta.vi[1]].pos[axis] +
                         mesh_->vertices[ta.vi[2]].pos[axis]) / 3.0;
            double cb = (mesh_->vertices[tb.vi[0]].pos[axis] +
                         mesh_->vertices[tb.vi[1]].pos[axis] +
                         mesh_->vertices[tb.vi[2]].pos[axis]) / 3.0;
            return ca < cb;
        });

    const int mid = (begin + end) / 2;

    // Recurse — note: buildNode may push_back into nodes_, so we must NOT hold
    // a reference/pointer to nodes_[node_idx] across the recursive calls.
    int left_child  = buildNode(face_ids, begin, mid);
    int right_child = buildNode(face_ids, mid,   end);

    // Re-index after potential reallocation.
    nodes_[node_idx].left  = left_child;
    nodes_[node_idx].right = right_child;

    return node_idx;
}

std::vector<uint32_t> BVH::queryZ(double z, double eps) const {
    std::vector<uint32_t> result;
    if (nodes_.empty()) return result;

    const double z_lo = z - eps;
    const double z_hi = z + eps;

    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);  // root is always node 0

    while (!stack.empty()) {
        const int idx = stack.back();
        stack.pop_back();

        const Node& node = nodes_[idx];

        if (!node.bounds.intersectsZSlab(z_lo, z_hi)) continue;

        if (node.left == -1) {
            // Leaf: collect face indices.
            const uint32_t end = node.face_begin + node.face_count;
            for (uint32_t i = node.face_begin; i < end; ++i) {
                result.push_back(leaf_faces_[i]);
            }
        } else {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// nearestTriangle
// ---------------------------------------------------------------------------

// Squared distance from point p to the axis-aligned bounding box.
static double aabbMinDistSq(const AABB& box, const Vec3d& p) {
    double d = 0.0;
    for (int k = 0; k < 3; ++k) {
        if      (p[k] < box.lo[k]) { double t = box.lo[k] - p[k]; d += t*t; }
        else if (p[k] > box.hi[k]) { double t = p[k] - box.hi[k]; d += t*t; }
    }
    return d;
}

// Squared distance from point p to triangle (A, B, C) — Eberly's method.
static double ptTriDistSq(const Vec3d& p,
                           const Vec3d& A, const Vec3d& B, const Vec3d& C)
{
    auto sub3 = [](const Vec3d& a, const Vec3d& b) -> Vec3d {
        return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
    };
    auto dot3 = [](const Vec3d& a, const Vec3d& b) {
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    };
    auto norm3sq = [&](const Vec3d& a){ return dot3(a,a); };

    Vec3d AB=sub3(B,A), AC=sub3(C,A), AP=sub3(p,A);
    double d1=dot3(AB,AP), d2=dot3(AC,AP);
    if(d1<=0&&d2<=0) return norm3sq(AP);

    Vec3d BP=sub3(p,B);
    double d3=dot3(AB,BP), d4=dot3(AC,BP);
    if(d3>=0&&d4<=d3) return norm3sq(BP);

    double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0){
        double v=d1/(d1-d3);
        Vec3d q={A[0]+v*AB[0],A[1]+v*AB[1],A[2]+v*AB[2]};
        return norm3sq(sub3(p,q));
    }

    Vec3d CP=sub3(p,C);
    double d5=dot3(AB,CP), d6=dot3(AC,CP);
    if(d6>=0&&d5<=d6) return norm3sq(CP);

    double vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0){
        double w=d2/(d2-d6);
        Vec3d q={A[0]+w*AC[0],A[1]+w*AC[1],A[2]+w*AC[2]};
        return norm3sq(sub3(p,q));
    }

    double va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){
        double w=(d4-d3)/((d4-d3)+(d5-d6));
        Vec3d BC=sub3(C,B);
        Vec3d q={B[0]+w*BC[0],B[1]+w*BC[1],B[2]+w*BC[2]};
        return norm3sq(sub3(p,q));
    }

    double denom=1.0/(va+vb+vc);
    double v=vb*denom, w=vc*denom;
    Vec3d q={A[0]+v*AB[0]+w*AC[0],A[1]+v*AB[1]+w*AC[1],A[2]+v*AB[2]+w*AC[2]};
    return norm3sq(sub3(p,q));
}

uint32_t BVH::nearestTriangle(const Vec3d& p, double& dist_out) const {
    if (nodes_.empty() || !mesh_) {
        dist_out = std::numeric_limits<double>::infinity();
        return 0;
    }

    double best_dist_sq = std::numeric_limits<double>::infinity();
    uint32_t best_face = 0;

    // Stack-based BVH traversal, pruning when AABB min-dist >= best.
    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(0);

    while (!stack.empty()) {
        const int idx = stack.back();
        stack.pop_back();

        const Node& node = nodes_[idx];
        if (aabbMinDistSq(node.bounds, p) >= best_dist_sq) continue;

        if (node.left == -1) {
            // Leaf: check all triangles.
            const uint32_t end = node.face_begin + node.face_count;
            for (uint32_t i = node.face_begin; i < end; ++i) {
                uint32_t fi = leaf_faces_[i];
                const auto& tri = mesh_->faces[fi];
                const Vec3d& A = mesh_->vertices[tri.vi[0]].pos;
                const Vec3d& B = mesh_->vertices[tri.vi[1]].pos;
                const Vec3d& C = mesh_->vertices[tri.vi[2]].pos;
                double dsq = ptTriDistSq(p, A, B, C);
                if (dsq < best_dist_sq) {
                    best_dist_sq = dsq;
                    best_face = fi;
                }
            }
        } else {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }

    dist_out = std::sqrt(best_dist_sq);
    return best_face;
}

} // namespace cfd
