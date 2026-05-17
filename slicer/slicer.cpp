#include "slicer.hpp"
#include "../geometry/bvh.hpp"
#include <cmath>
#include <fstream>
#include <unordered_map>
#include <cstdint>
#include <limits>

namespace cfd {

namespace {

inline int classify(double vz, double z) {
    // vertices exactly on the plane are treated as above to avoid duplicate segments
    return (vz >= z) ? 1 : -1;
}

inline Vec2d lerp_edge(const Vec3d& a, const Vec3d& b, double z) {
    double t = (z - a[2]) / (b[2] - a[2]);
    return { a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]) };
}

uint64_t snap_key(double x, double y) {
    int64_t ix = llround(x * 1e9);
    int64_t iy = llround(y * 1e9);
    // FNV-1a mix of both 64-bit integers
    uint64_t h = 14695981039346656037ULL;
    auto mix = [&](int64_t v) {
        for (int k = 0; k < 8; ++k) {
            h ^= static_cast<uint64_t>(v & 0xFF);
            h *= 1099511628211ULL;
            v >>= 8;
        }
    };
    mix(ix);
    mix(iy);
    return h;
}

using Segment = std::array<Vec2d, 2>;

std::vector<Segment> intersect_mesh(const SurfaceManifold& m, double z,
                                    const BVH& bvh) {
    std::vector<Segment> segs;

    auto candidate_faces = bvh.queryZ(z);
    segs.reserve(candidate_faces.size() / 2);

    for (uint32_t fi : candidate_faces) {
        const auto& tri = m.faces[fi];
        const Vec3d& p0 = m.vertices[tri.vi[0]].pos;
        const Vec3d& p1 = m.vertices[tri.vi[1]].pos;
        const Vec3d& p2 = m.vertices[tri.vi[2]].pos;

        int c0 = classify(p0[2], z);
        int c1 = classify(p1[2], z);
        int c2 = classify(p2[2], z);

        if (c0 == c1 && c1 == c2) continue;

        // collect the (exactly 2) crossing points from the three edges
        Vec2d pts[2];
        int n = 0;

        auto check_edge = [&](const Vec3d& a, const Vec3d& b, int ca, int cb) {
            if (ca == cb) return;
            pts[n++] = lerp_edge(a, b, z);
        };

        check_edge(p0, p1, c0, c1);
        check_edge(p1, p2, c1, c2);
        check_edge(p2, p0, c2, c0);

        // degenerate — shouldn't happen given classify treats on-plane as above
        if (n != 2) continue;

        segs.push_back({ pts[0], pts[1] });
    }
    return segs;
}

std::vector<Contour> chain_segments(const std::vector<Segment>& segs, double z) {
    struct HalfEdge { Vec2d from; Vec2d to; };
    std::unordered_map<uint64_t, HalfEdge> edge_map;
    edge_map.reserve(segs.size() * 2);

    for (const auto& s : segs) {
        uint64_t k0 = snap_key(s[0][0], s[0][1]);
        uint64_t k1 = snap_key(s[1][0], s[1][1]);
        edge_map[k0] = { s[0], s[1] };
        edge_map[k1] = { s[1], s[0] };
    }

    std::unordered_map<uint64_t, bool> used;
    used.reserve(edge_map.size());

    std::vector<Contour> contours;

    for (const auto& [start_key, he] : edge_map) {
        if (used[start_key]) continue;

        Contour c;
        c.z = z;

        uint64_t key = start_key;
        for (;;) {
            if (used[key]) break;
            used[key] = true;

            auto it = edge_map.find(key);
            if (it == edge_map.end()) break;

            const Vec2d& pt = it->second.from;

            // skip duplicate consecutive points
            if (c.points.empty() ||
                c.points.back()[0] != pt[0] || c.points.back()[1] != pt[1]) {
                c.points.push_back(pt);
            }

            key = snap_key(it->second.to[0], it->second.to[1]);
        }

        if (c.points.size() >= 2) {
            // remove closing duplicate if present
            if (c.points.front()[0] == c.points.back()[0] &&
                c.points.front()[1] == c.points.back()[1]) {
                c.points.pop_back();
            }
            if (c.points.size() >= 2)
                contours.push_back(std::move(c));
        }
    }

    return contours;
}

} // anonymous namespace

SliceResult Slicer::sliceAtZ(const SurfaceManifold& m, double z) const {
    BVH bvh;
    bvh.build(m);
    auto segs = intersect_mesh(m, z, bvh);
    SliceResult result;
    result.z_level = z;
    result.contours = chain_segments(segs, z);
    return result;
}

std::vector<SliceResult> Slicer::sliceRange(const SurfaceManifold& m,
                                              double z_min, double z_max,
                                              double z_interval) const {
    if (z_min == 0.0 && z_max == 0.0) {
        z_min = std::numeric_limits<double>::max();
        z_max = std::numeric_limits<double>::lowest();
        for (const auto& v : m.vertices) {
            if (v.pos[2] < z_min) z_min = v.pos[2];
            if (v.pos[2] > z_max) z_max = v.pos[2];
        }
    }

    std::vector<SliceResult> results;
    if (z_interval <= 0.0 || z_min >= z_max) return results;

    for (double z = z_min; z <= z_max + z_interval * 1e-9; z += z_interval) {
        if (z > z_max) z = z_max;
        results.push_back(sliceAtZ(m, z));
        if (z == z_max) break;
    }
    return results;
}

void Slicer::exportSVG(const std::vector<SliceResult>& slices,
                        const std::string& path) const {
    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::lowest();
    double ymax = std::numeric_limits<double>::lowest();

    for (const auto& sr : slices) {
        for (const auto& c : sr.contours) {
            for (const auto& p : c.points) {
                if (p[0] < xmin) xmin = p[0];
                if (p[1] < ymin) ymin = p[1];
                if (p[0] > xmax) xmax = p[0];
                if (p[1] > ymax) ymax = p[1];
            }
        }
    }

    if (xmin > xmax) { xmin = 0; ymin = 0; xmax = 1; ymax = 1; }

    double pad = (xmax - xmin + ymax - ymin) * 0.02;
    double vx = xmin - pad;
    double vy = ymin - pad;
    double vw = (xmax - xmin) + 2 * pad;
    double vh = (ymax - ymin) + 2 * pad;

    std::ofstream f(path);
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
    f << "width=\"800\" height=\"800\" ";
    f << "viewBox=\"" << vx << " " << vy << " " << vw << " " << vh << "\">\n";

    for (const auto& sr : slices) {
        for (const auto& c : sr.contours) {
            if (c.points.empty()) continue;
            f << "  <polyline fill=\"none\" stroke=\"black\" stroke-width=\""
              << vw * 0.001 << "\" points=\"";
            for (const auto& p : c.points) {
                f << p[0] << "," << p[1] << " ";
            }
            // close the loop
            f << c.points[0][0] << "," << c.points[0][1];
            f << "\"/>\n";
        }
    }

    f << "</svg>\n";
}

void Slicer::exportPolylines(const std::vector<SliceResult>& slices,
                               const std::string& path) const {
    std::ofstream f(path);
    for (const auto& sr : slices) {
        f << "# Z=" << sr.z_level << "\n";
        for (const auto& c : sr.contours) {
            for (const auto& p : c.points) {
                f << p[0] << " " << p[1] << "\n";
            }
            f << "\n";
        }
    }
}

} // namespace cfd
