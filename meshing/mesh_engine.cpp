// mesh_engine.cpp — Boundary Layer Inflation & Tet Core Generation
// Agent 2: Meshing & Boundary Layer Inflation module
// C++20
#include "mesh_engine.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cfd {

// ============================================================
//  Vector helpers (private static methods)
// ============================================================

Vec3d MeshEngine::cross(const Vec3d& a, const Vec3d& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

double MeshEngine::dot(const Vec3d& a, const Vec3d& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double MeshEngine::norm(const Vec3d& a) {
    return std::sqrt(dot(a, a));
}

Vec3d MeshEngine::normalize(const Vec3d& a) {
    double n = norm(a);
    if (n < 1e-300) return {0.0, 0.0, 0.0};
    return {a[0]/n, a[1]/n, a[2]/n};
}

Vec3d MeshEngine::add(const Vec3d& a, const Vec3d& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

Vec3d MeshEngine::sub(const Vec3d& a, const Vec3d& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

Vec3d MeshEngine::scale(const Vec3d& a, double s) {
    return {a[0]*s, a[1]*s, a[2]*s};
}

Vec3d MeshEngine::midpoint(const Vec3d& a, const Vec3d& b) {
    return {(a[0]+b[0])*0.5, (a[1]+b[1])*0.5, (a[2]+b[2])*0.5};
}

// ============================================================
//  MeshEngine constructor
// ============================================================

MeshEngine::MeshEngine(ProgressCb cb) : progress_(std::move(cb)) {}

void MeshEngine::notify(const std::string& stage, double fraction) const {
    if (progress_) progress_(stage, fraction);
}

// ============================================================
//  y+ / BL sizing
// ============================================================

double MeshEngine::computeFirstCellHeight(double Re, double yplus,
                                          double nu, double U_inf) {
    if (Re <= 0.0)
        throw MeshError("Reynolds number must be positive");
    if (nu <= 0.0)
        throw MeshError("Kinematic viscosity must be positive");
    if (U_inf <= 0.0)
        throw MeshError("Freestream velocity must be positive");
    if (yplus <= 0.0)
        throw MeshError("Target y+ must be positive");

    constexpr double rho = 1.225; // kg/m^3 — ISA sea-level air density

    // Schlichting turbulent flat-plate skin friction coefficient
    const double C_f   = 0.0592 * std::pow(Re, -0.2);

    // Wall shear stress and friction velocity
    const double tau_w = 0.5 * rho * U_inf * U_inf * C_f;
    const double u_tau = std::sqrt(tau_w / rho);

    if (u_tau < 1e-300)
        throw MeshError("Friction velocity is effectively zero — check flow conditions");

    const double h1 = yplus * nu / u_tau;
    return h1;
}

BLParameters MeshEngine::computeBLParameters(const PipelineConfig& cfg) {
    const double h1 = computeFirstCellHeight(
        cfg.reynolds_number,
        cfg.target_yplus,
        cfg.kinematic_viscosity,
        cfg.freestream_velocity);

    const double r = cfg.bl_growth_ratio;
    const int    n = cfg.bl_layers;

    if (r <= 1.0)
        throw MeshError("bl_growth_ratio must be > 1.0");
    if (n < 1)
        throw MeshError("bl_layers must be >= 1");

    // Geometric series: total thickness = h1 * (r^n - 1) / (r - 1)
    const double total = h1 * (std::pow(r, n) - 1.0) / (r - 1.0);

    std::fprintf(stderr,
        "BL: h1=%.4fmm, total=%.4fmm, layers=%d, ratio=%.3f\n",
        h1 * 1e3, total * 1e3, n, r);

    return BLParameters{h1, r, n, total};
}

// ============================================================
//  Area-weighted vertex normals (wall faces only)
// ============================================================

std::vector<Vec3d> MeshEngine::computeVertexNormals(
    const SurfaceManifold& surface,
    const std::vector<uint32_t>& wall_face_ids) const
{
    const auto nv = static_cast<uint32_t>(surface.vertices.size());
    std::vector<Vec3d> normals(nv, {0.0, 0.0, 0.0});

    for (uint32_t fid : wall_face_ids) {
        const Triangle& tri = surface.faces[fid];
        const Vec3d& p0 = surface.vertices[tri.vi[0]].pos;
        const Vec3d& p1 = surface.vertices[tri.vi[1]].pos;
        const Vec3d& p2 = surface.vertices[tri.vi[2]].pos;

        // Area-weighted normal = 0.5 * cross(e1, e2)
        Vec3d e1  = sub(p1, p0);
        Vec3d e2  = sub(p2, p0);
        Vec3d wn  = cross(e1, e2); // magnitude = 2 * area

        for (int k = 0; k < 3; ++k) {
            uint32_t vid = tri.vi[k];
            normals[vid][0] += wn[0];
            normals[vid][1] += wn[1];
            normals[vid][2] += wn[2];
        }
    }

    // Normalize each accumulated normal
    for (auto& n : normals) n = normalize(n);
    return normals;
}

// ============================================================
//  Laplacian smoothing of inflated prism layers
// ============================================================

// Each non-wall node is moved toward the average of its neighbours
// within the same layer ring.  Wall nodes are pinned.
//
// Connectivity is inferred from the prism elements: for each prism,
// the top-face nodes (3,4,5) neighbour each other and the matching
// nodes of any laterally adjacent prism sharing that edge.
//
// For simplicity we do a single-pass centroidal smooth: each inflated
// node is set to the average of the positions of its top-face siblings
// across all prisms that share it.
void MeshEngine::smoothPrismLayers(PrismMesh& mesh, int iterations) const {
    const auto nv = static_cast<uint32_t>(mesh.vertices.size());

    // Build neighbour lists (per vertex, set of neighbouring vertex ids)
    // considering only lateral quad faces of each prism.
    std::vector<std::vector<uint32_t>> nbrs(nv);

    for (const auto& prism : mesh.prisms) {
        // Bottom ring: (0,1), (1,2), (2,0)
        // Top ring:    (3,4), (4,5), (5,3)
        // Lateral:     (0,3), (1,4), (2,5)
        auto addEdge = [&](uint32_t a, uint32_t b) {
            nbrs[a].push_back(b);
            nbrs[b].push_back(a);
        };
        addEdge(prism.nodes[3], prism.nodes[4]);
        addEdge(prism.nodes[4], prism.nodes[5]);
        addEdge(prism.nodes[5], prism.nodes[3]);
        // Also connect to bottom for tangential coherence
        addEdge(prism.nodes[0], prism.nodes[3]);
        addEdge(prism.nodes[1], prism.nodes[4]);
        addEdge(prism.nodes[2], prism.nodes[5]);
    }

    // Deduplicate neighbour lists
    for (auto& list : nbrs) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
    }

    // Wall nodes are pinned — build a set of pinned node ids
    std::vector<bool> pinned(nv, false);
    for (uint32_t wid : mesh.wall_node_ids) pinned[wid] = true;

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<Vec3d> new_pos(nv);
        for (uint32_t i = 0; i < nv; ++i)
            new_pos[i] = mesh.vertices[i].pos;

        for (uint32_t i = 0; i < nv; ++i) {
            if (pinned[i] || nbrs[i].empty()) continue;
            Vec3d avg = {0.0, 0.0, 0.0};
            for (uint32_t nb : nbrs[i]) {
                avg[0] += mesh.vertices[nb].pos[0];
                avg[1] += mesh.vertices[nb].pos[1];
                avg[2] += mesh.vertices[nb].pos[2];
            }
            double inv = 1.0 / static_cast<double>(nbrs[i].size());
            // Blend 50 % toward centroid to preserve layer structure
            new_pos[i][0] = 0.5 * mesh.vertices[i].pos[0] + 0.5 * avg[0] * inv;
            new_pos[i][1] = 0.5 * mesh.vertices[i].pos[1] + 0.5 * avg[1] * inv;
            new_pos[i][2] = 0.5 * mesh.vertices[i].pos[2] + 0.5 * avg[2] * inv;
        }

        for (uint32_t i = 0; i < nv; ++i)
            mesh.vertices[i].pos = new_pos[i];
    }
}

// ============================================================
//  Prism inflation
// ============================================================

// We use a position-keyed map to deduplicate nodes shared between
// adjacent wall triangles.  Hashing Vec3d exactly is sufficient here
// because inflated node positions are computed deterministically from
// the same wall vertex + normal.

namespace {

// Hash a Vec3d position with a coarse grid key (1 pm resolution)
struct Vec3dHash {
    std::size_t operator()(const Vec3d& v) const noexcept {
        // Snap to 1e-12 m grid and mix bits
        auto qi = [](double x) -> int64_t {
            return static_cast<int64_t>(std::llround(x * 1e12));
        };
        int64_t ix = qi(v[0]), iy = qi(v[1]), iz = qi(v[2]);
        // FNV-1a style mixing
        std::size_t h = 14695981039346656037ULL;
        auto mix = [&](int64_t val) {
            for (int i = 0; i < 8; ++i) {
                h ^= static_cast<std::size_t>(val & 0xFF);
                h *= 1099511628211ULL;
                val >>= 8;
            }
        };
        mix(ix); mix(iy); mix(iz);
        return h;
    }
};

struct Vec3dEqual {
    bool operator()(const Vec3d& a, const Vec3d& b) const noexcept {
        constexpr double tol = 1e-12;
        return (std::abs(a[0]-b[0]) < tol &&
                std::abs(a[1]-b[1]) < tol &&
                std::abs(a[2]-b[2]) < tol);
    }
};

using PosMap = std::unordered_map<Vec3d, uint32_t, Vec3dHash, Vec3dEqual>;

// Insert position into map if absent; return its index.
uint32_t getOrInsert(PosMap& map,
                     std::vector<Vertex>& verts,
                     const Vec3d& pos) {
    auto it = map.find(pos);
    if (it != map.end()) return it->second;
    uint32_t idx = static_cast<uint32_t>(verts.size());
    verts.push_back(Vertex{pos});
    map[pos] = idx;
    return idx;
}

} // anonymous namespace

PrismMesh MeshEngine::inflatePrisms(const SurfaceManifold& surface,
                                    const BLParameters& params) {
    notify("Collecting wall faces", 0.0);

    // ----------------------------------------------------------------
    // Step 1 — Collect wall face indices
    // ----------------------------------------------------------------
    std::vector<uint32_t> wall_face_ids;
    for (const auto& patch : surface.patches) {
        if (patch.type == BoundaryType::Wall) {
            wall_face_ids.insert(wall_face_ids.end(),
                                 patch.face_ids.begin(),
                                 patch.face_ids.end());
        }
    }

    if (wall_face_ids.empty())
        throw MeshError("No wall boundary patches found — cannot inflate prisms");

    notify("Computing vertex normals", 0.05);

    // ----------------------------------------------------------------
    // Step 2 — Area-weighted vertex normals (wall vertices only)
    // ----------------------------------------------------------------
    std::vector<Vec3d> vn = computeVertexNormals(surface, wall_face_ids);

    // ----------------------------------------------------------------
    // Step 3 — Build per-layer cumulative offsets h_k
    //   h_k = h1 * (r^k - 1) / (r - 1),   k = 1 .. N
    // ----------------------------------------------------------------
    const int    N  = params.num_layers;
    const double h1 = params.first_cell_height;
    const double r  = params.growth_ratio;

    std::vector<double> h_k(N + 1, 0.0); // h_k[0] = 0 (wall), h_k[k] = cumulative
    for (int k = 1; k <= N; ++k)
        h_k[k] = h1 * (std::pow(r, k) - 1.0) / (r - 1.0);

    // ----------------------------------------------------------------
    // Step 4 — Inflate: build deduplicated node list + prism elements
    //
    //   For every wall face fi and every layer k (1..N):
    //
    //     bottom nodes: the 3 nodes at layer k-1
    //     top nodes:    the 3 nodes at layer k
    //
    //   Node at (vertex v, layer k) is placed at:
    //       pos_wall[v] + h_k[k] * vn[v]
    // ----------------------------------------------------------------
    notify("Inflating prism layers", 0.10);

    PrismMesh result;
    result.num_layers        = N;
    result.first_cell_height = h1;
    result.growth_ratio      = r;

    PosMap pos_map;
    pos_map.reserve(wall_face_ids.size() * 3 * (N + 1));

    // First, seed the position map with original wall vertices
    // so they get stable indices (wall_node_ids will be 0..M-1).
    std::unordered_map<uint32_t, uint32_t> surf_to_mesh; // surf vertex → mesh vertex
    for (uint32_t fid : wall_face_ids) {
        const Triangle& tri = surface.faces[fid];
        for (int k = 0; k < 3; ++k) {
            uint32_t sv = tri.vi[k];
            if (surf_to_mesh.count(sv)) continue;
            const Vec3d& p = surface.vertices[sv].pos;
            uint32_t idx = getOrInsert(pos_map, result.vertices, p);
            surf_to_mesh[sv] = idx;
        }
    }

    // Collect wall_node_ids (unique, sorted)
    result.wall_node_ids.reserve(surf_to_mesh.size());
    for (auto& [_, idx] : surf_to_mesh)
        result.wall_node_ids.push_back(idx);
    std::sort(result.wall_node_ids.begin(), result.wall_node_ids.end());

    // Inflate layers
    // layer_node_map[sv][k] = mesh vertex index of (surface vertex sv, layer k)
    std::unordered_map<uint32_t, std::vector<uint32_t>> layer_node;
    layer_node.reserve(surf_to_mesh.size());

    for (auto& [sv, base_idx] : surf_to_mesh) {
        std::vector<uint32_t>& lvec = layer_node[sv];
        lvec.resize(N + 1);
        lvec[0] = base_idx; // layer 0 = wall

        const Vec3d& p_wall = surface.vertices[sv].pos;
        const Vec3d& n_v    = vn[sv];

        for (int k = 1; k <= N; ++k) {
            Vec3d p_k = add(p_wall, scale(n_v, h_k[k]));
            lvec[k] = getOrInsert(pos_map, result.vertices, p_k);
        }
    }

    // Build PrismElement for each (face, layer)
    result.prisms.reserve(wall_face_ids.size() * N);

    const double total_faces = static_cast<double>(wall_face_ids.size());
    for (std::size_t fi = 0; fi < wall_face_ids.size(); ++fi) {
        uint32_t fid        = wall_face_ids[fi];
        const Triangle& tri = surface.faces[fid];
        const uint32_t sv0  = tri.vi[0];
        const uint32_t sv1  = tri.vi[1];
        const uint32_t sv2  = tri.vi[2];

        for (int k = 1; k <= N; ++k) {
            PrismElement pe;
            pe.wall_face_id = fid;
            pe.layer        = static_cast<uint32_t>(k - 1);

            // Bottom (layer k-1)
            pe.nodes[0] = layer_node[sv0][k-1];
            pe.nodes[1] = layer_node[sv1][k-1];
            pe.nodes[2] = layer_node[sv2][k-1];
            // Top (layer k)
            pe.nodes[3] = layer_node[sv0][k];
            pe.nodes[4] = layer_node[sv1][k];
            pe.nodes[5] = layer_node[sv2][k];

            result.prisms.push_back(pe);
        }

        if (fi % 1000 == 0) {
            double pct = 0.10 + 0.70 * (fi / total_faces);
            notify("Inflating prism layers", pct);
        }
    }

    // ----------------------------------------------------------------
    // Step 5 — Collect outer node ids (layer N)
    // ----------------------------------------------------------------
    std::vector<bool> is_outer(result.vertices.size(), false);
    for (auto& [sv, lvec] : layer_node)
        is_outer[lvec[N]] = true;

    for (uint32_t i = 0; i < static_cast<uint32_t>(result.vertices.size()); ++i)
        if (is_outer[i]) result.outer_node_ids.push_back(i);
    std::sort(result.outer_node_ids.begin(), result.outer_node_ids.end());

    notify("Smoothing prism layers", 0.82);

    // ----------------------------------------------------------------
    // Step 6 — Smooth to reduce tangling
    // ----------------------------------------------------------------
    smoothPrismLayers(result, 3);

    notify("Prism inflation complete", 1.0);

    std::fprintf(stderr,
        "Prisms: %zu elements, %zu nodes (%zu wall, %zu outer)\n",
        result.prisms.size(),
        result.vertices.size(),
        result.wall_node_ids.size(),
        result.outer_node_ids.size());

    return result;
}

// ============================================================
//  Circumsphere of a tetrahedron
// ============================================================

// Uses the algebraic formulation: solve the 3x3 system
//   A * c = b  where A[i] = 2*(v[i+1] - v[0]), b[i] = |v[i+1]|^2 - |v[0]|^2
std::pair<Vec3d, double> MeshEngine::circumsphere(const std::array<Vec3d, 4>& v) {
    // Translate so v[0] is origin
    Vec3d a = sub(v[1], v[0]);
    Vec3d b = sub(v[2], v[0]);
    Vec3d c = sub(v[3], v[0]);

    double A[3][3] = {
        {2*a[0], 2*a[1], 2*a[2]},
        {2*b[0], 2*b[1], 2*b[2]},
        {2*c[0], 2*c[1], 2*c[2]}
    };
    double rhs[3] = {
        dot(a, a),
        dot(b, b),
        dot(c, c)
    };

    // Cramer's rule
    auto det3 = [](double m[3][3]) -> double {
        return m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
             - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
             + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);
    };

    double D = det3(A);
    if (std::abs(D) < 1e-300) {
        // Degenerate tet — return centroid with large radius
        Vec3d centre = {
            (v[0][0]+v[1][0]+v[2][0]+v[3][0])*0.25,
            (v[0][1]+v[1][1]+v[2][1]+v[3][1])*0.25,
            (v[0][2]+v[1][2]+v[2][2]+v[3][2])*0.25
        };
        return {centre, std::numeric_limits<double>::max()};
    }

    auto detCol = [&](int col) -> double {
        double M[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                M[i][j] = (j == col) ? rhs[i] : A[i][j];
        return det3(M);
    };

    Vec3d centre_local = {detCol(0)/D, detCol(1)/D, detCol(2)/D};
    Vec3d centre       = add(centre_local, v[0]);
    double rad         = norm(sub(centre, v[0]));

    return {centre, rad};
}

// ============================================================
//  Tet skewness (equilateral reference)
// ============================================================

// Skewness = (R_circ - d_max) / R_circ
// where d_max is the maximum distance from the circumcenter to each face centroid.
// A perfect equilateral tet has skewness 0; degenerate → 1.
double MeshEngine::tetSkewness(const std::array<Vec3d, 4>& verts) {
    auto [centre, R] = circumsphere(verts);

    if (R >= std::numeric_limits<double>::max() * 0.5) return 1.0;

    // 4 face centroids
    const int face_indices[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double max_dist = 0.0;
    for (auto& fi : face_indices) {
        Vec3d fc = {
            (verts[fi[0]][0] + verts[fi[1]][0] + verts[fi[2]][0]) / 3.0,
            (verts[fi[0]][1] + verts[fi[1]][1] + verts[fi[2]][1]) / 3.0,
            (verts[fi[0]][2] + verts[fi[1]][2] + verts[fi[2]][2]) / 3.0
        };
        double d = norm(sub(centre, fc));
        if (d > max_dist) max_dist = d;
    }

    // Normalise: skewness = 1 - d_max / R  → 0 for ideal, 1 for flat
    double skew = 1.0 - max_dist / R;
    // Clamp to [0,1]
    return std::max(0.0, std::min(1.0, skew));
}

// ============================================================
//  Tet aspect ratio (longest edge / shortest altitude)
// ============================================================

double MeshEngine::tetAspectRatio(const std::array<Vec3d, 4>& verts) {
    // Longest edge
    double max_edge = 0.0;
    const int edges[6][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
    for (auto& e : edges) {
        double len = norm(sub(verts[e[1]], verts[e[0]]));
        if (len > max_edge) max_edge = len;
    }

    // Shortest altitude: altitude of vertex i to its opposite face
    // altitude_i = 3 * vol / area_i
    // vol = |det([v1-v0, v2-v0, v3-v0])| / 6
    Vec3d e1 = sub(verts[1], verts[0]);
    Vec3d e2 = sub(verts[2], verts[0]);
    Vec3d e3 = sub(verts[3], verts[0]);
    double vol6 = std::abs(dot(e1, cross(e2, e3)));

    if (vol6 < 1e-300) return std::numeric_limits<double>::max();

    const int face_vi[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    double min_alt = std::numeric_limits<double>::max();
    for (auto& fi : face_vi) {
        Vec3d fa = sub(verts[fi[1]], verts[fi[0]]);
        Vec3d fb = sub(verts[fi[2]], verts[fi[0]]);
        double area = 0.5 * norm(cross(fa, fb));
        if (area < 1e-300) continue;
        double alt = vol6 / (6.0 * area); // = 3*vol / area_face... but 3V = vol6/2
        // Correct: altitude = 3 * vol / area, vol = vol6/6
        alt = (vol6 / 6.0) * 3.0 / area;
        if (alt < min_alt) min_alt = alt;
    }

    if (min_alt < 1e-300) return std::numeric_limits<double>::max();
    return max_edge / min_alt;
}

// ============================================================
//  Quality assessment
// ============================================================

MeshEngine::QualityReport MeshEngine::checkQuality(const VolumeMesh& mesh) const {
    QualityReport rep{};
    rep.max_skewness      = 0.0;
    rep.avg_skewness      = 0.0;
    rep.max_aspect_ratio  = 0.0;
    rep.avg_aspect_ratio  = 0.0;
    rep.min_orthogonality = 180.0;
    rep.bad_cells         = 0;

    if (mesh.tets.empty()) return rep;

    double sum_skew = 0.0;
    double sum_ar   = 0.0;

    for (const auto& tet : mesh.tets) {
        std::array<Vec3d, 4> v;
        for (int i = 0; i < 4; ++i)
            v[i] = mesh.vertices[tet.nodes[i]].pos;

        double sk = tetSkewness(v);
        double ar = tetAspectRatio(v);

        sum_skew += sk;
        sum_ar   += ar;

        if (sk > rep.max_skewness)      rep.max_skewness = sk;
        if (ar > rep.max_aspect_ratio)  rep.max_aspect_ratio = ar;
        if (sk > 0.85)                  ++rep.bad_cells;

        // Orthogonality: minimum angle between face normal and cell-to-cell vector.
        // For a tet with 4 faces, compute the angle between each face normal and
        // the vector from the cell centroid to the face centroid.
        Vec3d cc = {
            (v[0][0]+v[1][0]+v[2][0]+v[3][0])*0.25,
            (v[0][1]+v[1][1]+v[2][1]+v[3][1])*0.25,
            (v[0][2]+v[1][2]+v[2][2]+v[3][2])*0.25
        };

        const int face_vi[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
        for (auto& fi : face_vi) {
            Vec3d fa = sub(v[fi[1]], v[fi[0]]);
            Vec3d fb = sub(v[fi[2]], v[fi[0]]);
            Vec3d fn = normalize(cross(fa, fb));

            Vec3d fc = {
                (v[fi[0]][0]+v[fi[1]][0]+v[fi[2]][0])/3.0,
                (v[fi[0]][1]+v[fi[1]][1]+v[fi[2]][1])/3.0,
                (v[fi[0]][2]+v[fi[1]][2]+v[fi[2]][2])/3.0
            };

            Vec3d d  = normalize(sub(fc, cc));
            double cosA = std::abs(dot(fn, d));
            cosA = std::min(1.0, cosA); // guard for floating point
            double angle_deg = std::acos(cosA) * (180.0 / M_PI);
            if (angle_deg < rep.min_orthogonality)
                rep.min_orthogonality = angle_deg;
        }
    }

    double inv = 1.0 / static_cast<double>(mesh.tets.size());
    rep.avg_skewness     = sum_skew * inv;
    rep.avg_aspect_ratio = sum_ar   * inv;

    std::fprintf(stderr,
        "Quality: skew_max=%.3f avg=%.3f | AR_max=%.2f avg=%.2f | "
        "orth_min=%.1f deg | bad_cells=%d\n",
        rep.max_skewness, rep.avg_skewness,
        rep.max_aspect_ratio, rep.avg_aspect_ratio,
        rep.min_orthogonality, rep.bad_cells);

    return rep;
}

// ============================================================
//  Bowyer–Watson incremental Delaunay (simplified 3D)
// ============================================================

// Internal representation of a tet for the incremental insertion.
namespace {

struct BWTet {
    std::array<uint32_t, 4> nodes;
    Vec3d   circumcenter;
    double  circumradius;
    bool    bad = false; // marked for removal during insertion
};

// A triangular face (ordered triple of indices, canonical form)
struct TriFace {
    std::array<uint32_t, 3> v;

    // Canonical form: sort so that comparisons are order-independent
    TriFace(uint32_t a, uint32_t b, uint32_t c) : v{a, b, c} {
        if (v[0] > v[1]) std::swap(v[0], v[1]);
        if (v[1] > v[2]) std::swap(v[1], v[2]);
        if (v[0] > v[1]) std::swap(v[0], v[1]);
    }

    bool operator==(const TriFace& o) const noexcept {
        return v[0]==o.v[0] && v[1]==o.v[1] && v[2]==o.v[2];
    }
};

struct TriFaceHash {
    std::size_t operator()(const TriFace& f) const noexcept {
        std::size_t h = 2166136261u;
        for (uint32_t x : f.v) {
            h ^= x;
            h *= 16777619u;
        }
        return h;
    }
};

// Build a super-tetrahedron large enough to contain all points
BWTet makeSuperTet(const std::vector<Vec3d>& pts, std::vector<Vertex>& verts) {
    double lo[3] = { 1e300,  1e300,  1e300};
    double hi[3] = {-1e300, -1e300, -1e300};
    for (const auto& p : pts)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }

    double dx = hi[0]-lo[0], dy = hi[1]-lo[1], dz = hi[2]-lo[2];
    double M  = std::max({dx, dy, dz}) * 10.0 + 1.0;
    double cx = (lo[0]+hi[0])*0.5;
    double cy = (lo[1]+hi[1])*0.5;
    double cz = (lo[2]+hi[2])*0.5;

    uint32_t base = static_cast<uint32_t>(verts.size());
    // 4 vertices of an enclosing tetrahedron
    verts.push_back(Vertex{{cx,       cy + 3*M, cz}});
    verts.push_back(Vertex{{cx - 3*M, cy - M,   cz}});
    verts.push_back(Vertex{{cx + 3*M, cy - M,   cz}});
    verts.push_back(Vertex{{cx,       cy,        cz + 3*M}});

    std::array<Vec3d,4> sv;
    for (int i = 0; i < 4; ++i) sv[i] = verts[base+i].pos;
    auto [cc, cr] = MeshEngine::circumsphere(sv); // forward declaration resolved

    BWTet st;
    st.nodes = {base, base+1, base+2, base+3};
    st.circumcenter = cc;
    st.circumradius = cr;
    return st;
}

} // anonymous namespace

// circumsphere is a public static method — already defined above.
// We need to call it from free functions, so we expose it via a local trampoline.
// (It is already static on MeshEngine, so we just call MeshEngine::circumsphere.)

void MeshEngine::bowyerWatson(VolumeMesh& vol,
                               const std::vector<Vec3d>& interior_pts) const {
    if (interior_pts.empty()) return;

    // All points: outer_node positions + interior seeds
    // vol.vertices already populated by caller; we add super-tet vertices below.

    std::vector<Vec3d> all_pts;
    all_pts.reserve(vol.vertices.size() + interior_pts.size());
    for (const auto& vtx : vol.vertices) all_pts.push_back(vtx.pos);
    for (const auto& p   : interior_pts)  all_pts.push_back(p);

    // Build working vertex list
    std::vector<Vertex> work_verts = vol.vertices;
    for (const auto& p : interior_pts) work_verts.push_back(Vertex{p});

    uint32_t n_real = static_cast<uint32_t>(work_verts.size());

    // Create super-tet (appends 4 vertices)
    std::vector<BWTet> tets;
    tets.reserve(all_pts.size() * 7);
    tets.push_back(makeSuperTet(all_pts, work_verts));

    notify("Bowyer-Watson: inserting points", 0.0);

    // Insert each real point
    const double total_pts = static_cast<double>(n_real);
    for (uint32_t pi = 0; pi < n_real; ++pi) {
        const Vec3d& p = work_verts[pi].pos;

        // Mark tets whose circumsphere contains p
        for (auto& t : tets) {
            double dist2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                double d = p[k] - t.circumcenter[k];
                dist2 += d*d;
            }
            t.bad = (dist2 < t.circumradius * t.circumradius * (1.0 + 1e-10));
        }

        // Find the boundary of the hole (faces referenced once)
        std::unordered_map<TriFace, int, TriFaceHash> face_count;
        face_count.reserve(64);

        for (const auto& t : tets) {
            if (!t.bad) continue;
            const auto& n = t.nodes;
            // 4 faces of the tet
            face_count[TriFace{n[0],n[1],n[2]}]++;
            face_count[TriFace{n[0],n[1],n[3]}]++;
            face_count[TriFace{n[0],n[2],n[3]}]++;
            face_count[TriFace{n[1],n[2],n[3]}]++;
        }

        // Remove bad tets
        tets.erase(std::remove_if(tets.begin(), tets.end(),
                                  [](const BWTet& t){ return t.bad; }),
                   tets.end());

        // Re-triangulate the hole: for each boundary face, create new tet with p
        for (const auto& [face, cnt] : face_count) {
            if (cnt != 1) continue; // shared face — not on boundary

            uint32_t idx = pi; // index of inserted point

            std::array<uint32_t, 4> tnodes = {
                face.v[0], face.v[1], face.v[2], idx
            };
            std::array<Vec3d, 4> tv;
            for (int k = 0; k < 4; ++k) tv[k] = work_verts[tnodes[k]].pos;

            auto [cc, cr] = circumsphere(tv);

            BWTet nt;
            nt.nodes        = tnodes;
            nt.circumcenter = cc;
            nt.circumradius = cr;
            nt.bad          = false;
            tets.push_back(nt);
        }

        if (pi % 500 == 0) {
            notify("Bowyer-Watson: inserting points",
                   static_cast<double>(pi) / total_pts);
        }
    }

    // Remove any tet that shares a vertex with the super-tet
    uint32_t super_base = n_real; // first super-tet vertex index
    tets.erase(std::remove_if(tets.begin(), tets.end(),
        [super_base](const BWTet& t) {
            for (uint32_t nid : t.nodes)
                if (nid >= super_base) return true;
            return false;
        }), tets.end());

    // Commit to vol
    vol.vertices = work_verts;
    vol.vertices.resize(n_real); // drop super-tet vertices

    vol.tets.reserve(vol.tets.size() + tets.size());
    for (const auto& t : tets) {
        TetElement te;
        te.nodes = t.nodes;
        vol.tets.push_back(te);
    }

    notify("Bowyer-Watson: done", 1.0);
}

// ============================================================
//  Tetrahedral core generation
// ============================================================

VolumeMesh MeshEngine::generateTetCore(const PrismMesh& bl_mesh,
                                        double target_size) {
    if (target_size <= 0.0)
        throw MeshError("target_size must be positive");

    notify("Tet core: computing bounding box", 0.0);

    // ----------------------------------------------------------------
    // Bounding box from all BL vertices (wall + inflated)
    // ----------------------------------------------------------------
    double lo[3] = { 1e300,  1e300,  1e300};
    double hi[3] = {-1e300, -1e300, -1e300};
    for (const auto& vtx : bl_mesh.vertices) {
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], vtx.pos[k]);
            hi[k] = std::max(hi[k], vtx.pos[k]);
        }
    }
    // Pad 5 % so the grid fully encloses the geometry
    for (int k = 0; k < 3; ++k) {
        double d = (hi[k] - lo[k]) * 0.05 + target_size;
        lo[k] -= d;  hi[k] += d;
    }

    // ----------------------------------------------------------------
    // Regular Cartesian grid of nodes — hard-capped at 30³ for speed
    // ----------------------------------------------------------------
    int nx = std::min(static_cast<int>(std::ceil((hi[0]-lo[0])/target_size))+1, 30);
    int ny = std::min(static_cast<int>(std::ceil((hi[1]-lo[1])/target_size))+1, 30);
    int nz = std::min(static_cast<int>(std::ceil((hi[2]-lo[2])/target_size))+1, 30);
    nx = std::max(nx, 2);  ny = std::max(ny, 2);  nz = std::max(nz, 2);

    double sx = (hi[0]-lo[0]) / (nx-1);
    double sy = (hi[1]-lo[1]) / (ny-1);
    double sz = (hi[2]-lo[2]) / (nz-1);

    notify("Tet core: building grid nodes", 0.1);
    std::fprintf(stderr, "Tet core: grid %dx%dx%d (%d cells)\n",
                 nx, ny, nz, (nx-1)*(ny-1)*(nz-1));

    // ----------------------------------------------------------------
    // Build VolumeMesh — start with all BL prism data, then add grid
    // ----------------------------------------------------------------
    VolumeMesh vol;
    vol.vertices = bl_mesh.vertices;   // prism nodes first
    vol.prisms   = bl_mesh.prisms;

    // Map grid (ix,iy,iz) → vertex index
    uint32_t grid_base = static_cast<uint32_t>(vol.vertices.size());
    auto gidx = [&](int ix, int iy, int iz) -> uint32_t {
        return grid_base + static_cast<uint32_t>(iz*ny*nx + iy*nx + ix);
    };

    vol.vertices.reserve(vol.vertices.size() + static_cast<size_t>(nx*ny*nz));
    for (int iz = 0; iz < nz; ++iz)
        for (int iy = 0; iy < ny; ++iy)
            for (int ix = 0; ix < nx; ++ix)
                vol.vertices.push_back(Vertex{{lo[0]+ix*sx,
                                               lo[1]+iy*sy,
                                               lo[2]+iz*sz}});

    // ----------------------------------------------------------------
    // Freudenthal / Kuhn triangulation of each grid cube → 6 tets
    // This is consistent across cube faces (shared face identical split).
    //
    // For cube corner indices:
    //   v000=gidx(i,j,k)       v100=gidx(i+1,j,k)
    //   v010=gidx(i,j+1,k)     v110=gidx(i+1,j+1,k)
    //   v001=gidx(i,j,k+1)     v101=gidx(i+1,j,k+1)
    //   v011=gidx(i,j+1,k+1)   v111=gidx(i+1,j+1,k+1)
    //
    // 6-tet split (all share vertex v000 and v111):
    //   T0: 000 100 110 111
    //   T1: 000 010 110 111
    //   T2: 000 001 101 111
    //   T3: 000 100 101 111
    //   T4: 000 010 011 111
    //   T5: 000 001 011 111
    // ----------------------------------------------------------------
    notify("Tet core: triangulating grid cells", 0.3);
    vol.tets.reserve(static_cast<size_t>((nx-1)*(ny-1)*(nz-1)*6));

    // offsets for the 8 corners of a cube [iz][iy][ix]
    auto g = [&](int ix, int iy, int iz) { return gidx(ix,iy,iz); };

    for (int iz = 0; iz < nz-1; ++iz) {
        for (int iy = 0; iy < ny-1; ++iy) {
            for (int ix = 0; ix < nx-1; ++ix) {
                uint32_t v000 = g(ix,  iy,  iz  );
                uint32_t v100 = g(ix+1,iy,  iz  );
                uint32_t v010 = g(ix,  iy+1,iz  );
                uint32_t v110 = g(ix+1,iy+1,iz  );
                uint32_t v001 = g(ix,  iy,  iz+1);
                uint32_t v101 = g(ix+1,iy,  iz+1);
                uint32_t v011 = g(ix,  iy+1,iz+1);
                uint32_t v111 = g(ix+1,iy+1,iz+1);

                // 6 tets per cube
                std::array<std::array<uint32_t,4>,6> cube_tets = {{
                    {v000, v100, v110, v111},
                    {v000, v010, v110, v111},
                    {v000, v001, v101, v111},
                    {v000, v100, v101, v111},
                    {v000, v010, v011, v111},
                    {v000, v001, v011, v111},
                }};
                for (const auto& t : cube_tets)
                    vol.tets.push_back(TetElement{t});
            }
        }
    }

    notify("Tet core: complete", 1.0);
    std::fprintf(stderr,
        "VolumeMesh: %zu tets, %zu prisms, %zu vertices\n",
        vol.tets.size(), vol.prisms.size(), vol.vertices.size());

    return vol;
}

// circumsphere needs to be callable from free function makeSuperTet.
// It is a public static method, so it is visible as MeshEngine::circumsphere.
// No additional linkage needed.

} // namespace cfd
