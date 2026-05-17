#include "curvature.hpp"
#include <cmath>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cstdint>

namespace cfd {

namespace {

double dot3(const Vec3d& a, const Vec3d& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vec3d cross3(const Vec3d& a, const Vec3d& b) {
    return { a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0] };
}

double norm3(const Vec3d& a) {
    return std::sqrt(dot3(a, a));
}

Vec3d sub3(const Vec3d& a, const Vec3d& b) {
    return { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
}

// Edge key: encodes an undirected edge as a single uint64 with lo index in high bits.
uint64_t edgeKey(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

} // namespace

double CurvatureAnalyzer::cotWeight(const Vec3d& p, const Vec3d& q, const Vec3d& r) {
    const Vec3d u = sub3(q, p);
    const Vec3d v = sub3(r, p);
    const double d = dot3(u, v);
    const double c = norm3(cross3(u, v));
    // c == 0 means degenerate triangle; return 0 to skip its contribution.
    if (c == 0.0) return 0.0;
    return d / c;
}

std::array<uint8_t,3> CurvatureAnalyzer::jetColor(double t) {
    t = std::min(std::max(t, 0.0), 1.0);

    double r, g, b;
    if (t < 0.25) {
        const double s = t / 0.25;
        r = 0.0;
        g = s;
        b = 1.0;
    } else if (t < 0.5) {
        const double s = (t - 0.25) / 0.25;
        r = 0.0;
        g = 1.0;
        b = 1.0 - s;
    } else if (t < 0.75) {
        const double s = (t - 0.5) / 0.25;
        r = s;
        g = 1.0;
        b = 0.0;
    } else {
        const double s = (t - 0.75) / 0.25;
        r = 1.0;
        g = 1.0 - s;
        b = 0.0;
    }

    return {
        static_cast<uint8_t>(std::round(r * 255.0)),
        static_cast<uint8_t>(std::round(g * 255.0)),
        static_cast<uint8_t>(std::round(b * 255.0))
    };
}

CurvatureResult CurvatureAnalyzer::compute(const SurfaceManifold& m) const {
    const uint32_t nv = static_cast<uint32_t>(m.vertices.size());

    // --- adjacency build: vertex -> neighbour vertices (1-ring) ---
    std::unordered_map<uint32_t, std::vector<uint32_t>> adj;
    adj.reserve(nv);
    for (const Triangle& tri : m.faces) {
        for (int k = 0; k < 3; ++k) {
            const uint32_t a = tri.vi[k];
            const uint32_t b = tri.vi[(k+1)%3];
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
    }

    // --- per-edge data: both opposite vertices for the two triangles sharing it ---
    // Maps edgeKey -> pair of opposite vertex indices (filled as triangles are scanned).
    struct EdgeOpposites {
        uint32_t opp[2] = {UINT32_MAX, UINT32_MAX};
        int count = 0;
    };
    std::unordered_map<uint64_t, EdgeOpposites> edgeOpp;
    edgeOpp.reserve(m.faces.size() * 3);

    for (const Triangle& tri : m.faces) {
        for (int k = 0; k < 3; ++k) {
            const uint32_t vi = tri.vi[k];
            const uint32_t vj = tri.vi[(k+1)%3];
            const uint32_t vo = tri.vi[(k+2)%3]; // opposite vertex
            auto& eo = edgeOpp[edgeKey(vi, vj)];
            if (eo.count < 2) eo.opp[eo.count++] = vo;
        }
    }

    // --- cotangent-weight Laplace-Beltrami accumulation (single pass) ---
    std::vector<double> H(nv, 0.0);
    std::vector<double> A_mixed(nv, 0.0);
    std::vector<Vec3d> H_vec(nv, {0.0, 0.0, 0.0});

    for (const auto& [key, eo] : edgeOpp) {
        if (eo.count < 2) continue; // boundary edge — skip (no second triangle)

        const uint32_t i = static_cast<uint32_t>(key >> 32);
        const uint32_t j = static_cast<uint32_t>(key & 0xFFFFFFFF);

        const Vec3d& pi = m.vertices[i].pos;
        const Vec3d& pj = m.vertices[j].pos;
        const Vec3d& pa = m.vertices[eo.opp[0]].pos;
        const Vec3d& pb = m.vertices[eo.opp[1]].pos;

        const double cot_a = cotWeight(pa, pi, pj);
        const double cot_b = cotWeight(pb, pi, pj);
        const double w = cot_a + cot_b;

        const Vec3d diff_ij = sub3(pi, pj);
        const double edge_len_sq = dot3(diff_ij, diff_ij);

        // For mixed area: each edge contributes (1/8)*w*|e|^2 to both endpoint vertices.
        const double area_contrib = (1.0 / 8.0) * w * edge_len_sq;
        A_mixed[i] += area_contrib;
        A_mixed[j] += area_contrib;

        // Accumulate weighted displacement into H_vec for both endpoints.
        H_vec[i][0] += w * diff_ij[0];
        H_vec[i][1] += w * diff_ij[1];
        H_vec[i][2] += w * diff_ij[2];

        const Vec3d diff_ji = sub3(pj, pi);
        H_vec[j][0] += w * diff_ji[0];
        H_vec[j][1] += w * diff_ji[1];
        H_vec[j][2] += w * diff_ji[2];
    }

    for (uint32_t i = 0; i < nv; ++i) {
        const double A = std::max(A_mixed[i], 1e-15);
        H[i] = norm3(H_vec[i]) / (2.0 * A);
    }

    // --- 2-sigma clamp then normalise to [0,1] for jet colormap ---
    const double mean_H = std::accumulate(H.begin(), H.end(), 0.0) / static_cast<double>(nv);

    double variance = 0.0;
    for (const double h : H) {
        const double d = h - mean_H;
        variance += d * d;
    }
    variance /= static_cast<double>(nv);
    const double sigma = std::sqrt(variance);

    const double lo = mean_H - 2.0 * sigma;
    const double hi = mean_H + 2.0 * sigma;
    const double range = hi - lo;

    CurvatureResult result;
    result.mean_curvature = H;
    result.jet_colors.resize(nv);

    for (uint32_t i = 0; i < nv; ++i) {
        const double t = (range > 0.0) ? (H[i] - lo) / range : 0.5;
        result.jet_colors[i] = jetColor(t);
    }

    return result;
}

} // namespace cfd
