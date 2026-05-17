#include "geometry_engine.hpp"

#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <queue>
#include <cassert>
#include <cctype>
#include <stdexcept>

namespace cfd {

// ---------------------------------------------------------------------------
// Pimpl body (opaque state, reserved for future OCC / Gmsh handles)
// ---------------------------------------------------------------------------
struct GeometryEngine::Impl {
    // Placeholder for OpenCASCADE / Gmsh context handles.
    // Additional per-engine state can be added here without changing the ABI.
    bool gmsh_initialized = false;
};

// ---------------------------------------------------------------------------
// Constructor / destructor
// ---------------------------------------------------------------------------
GeometryEngine::GeometryEngine(ProgressCb cb)
    : impl_(std::make_unique<Impl>()), progress_(std::move(cb)) {}

GeometryEngine::~GeometryEngine() = default;

// ---------------------------------------------------------------------------
// Internal helper: fire progress callback if one was supplied
// ---------------------------------------------------------------------------
static void emit(const GeometryEngine::ProgressCb& cb,
                 const std::string& stage, double frac) {
    if (cb) cb(stage, frac);
}

// ---------------------------------------------------------------------------
// Spatial-hash helpers
// ---------------------------------------------------------------------------
static uint64_t hashCoord(double x, double y, double z, double inv_tol) {
    // Map each coordinate to an integer grid cell and pack into a 64-bit key.
    // Uses prime multipliers to reduce collisions.
    auto ix = static_cast<int64_t>(std::floor(x * inv_tol));
    auto iy = static_cast<int64_t>(std::floor(y * inv_tol));
    auto iz = static_cast<int64_t>(std::floor(z * inv_tol));

    // Fold negatives: shift into positive range with a large offset then mask.
    constexpr int64_t OFFSET = (1LL << 20);  // ±1 M cells per axis
    uint64_t ux = static_cast<uint64_t>(ix + OFFSET) & 0xFFFFFull;
    uint64_t uy = static_cast<uint64_t>(iy + OFFSET) & 0xFFFFFull;
    uint64_t uz = static_cast<uint64_t>(iz + OFFSET) & 0xFFFFFull;
    return (ux) | (uy << 20) | (uz << 40);
}

// ---------------------------------------------------------------------------
// mergeVertices
// ---------------------------------------------------------------------------
void GeometryEngine::mergeVertices(SurfaceManifold& m, double tol) {
    if (m.vertices.empty()) return;

    const double inv_tol = 1.0 / tol;
    std::unordered_map<uint64_t, uint32_t> grid;
    grid.reserve(m.vertices.size() * 2);

    // old-index → canonical new index
    std::vector<uint32_t> remap(m.vertices.size());
    std::vector<Vertex>   new_verts;
    new_verts.reserve(m.vertices.size());

    for (uint32_t i = 0; i < static_cast<uint32_t>(m.vertices.size()); ++i) {
        const auto& p = m.vertices[i].pos;
        uint64_t key  = hashCoord(p[0], p[1], p[2], inv_tol);

        if (auto it = grid.find(key); it != grid.end()) {
            remap[i] = it->second;
        } else {
            uint32_t new_idx = static_cast<uint32_t>(new_verts.size());
            new_verts.push_back(m.vertices[i]);
            grid[key] = new_idx;
            remap[i]  = new_idx;
        }
    }

    m.vertices = std::move(new_verts);

    // Remap triangle vertex indices
    for (auto& tri : m.faces) {
        for (auto& vi : tri.vi) {
            vi = remap[vi];
        }
    }
}

// ---------------------------------------------------------------------------
// fixWinding — BFS flood-fill from face 0
// ---------------------------------------------------------------------------
void GeometryEngine::fixWinding(SurfaceManifold& m) {
    if (m.faces.empty()) return;

    const uint32_t nf = static_cast<uint32_t>(m.faces.size());

    // Build half-edge adjacency: directed edge (a→b) → face index
    // Key: packed (a, b) with a < max(uint32_t)
    using EdgeKey = uint64_t;
    auto makeEdgeKey = [](uint32_t a, uint32_t b) -> EdgeKey {
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    std::unordered_map<EdgeKey, uint32_t> half_edge_face;
    half_edge_face.reserve(nf * 6);

    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& tri = m.faces[fi];
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            half_edge_face[makeEdgeKey(a, b)] = fi;
        }
    }

    std::vector<bool> visited(nf, false);
    std::vector<bool> flipped(nf, false);
    std::queue<uint32_t> bfs;

    bfs.push(0);
    visited[0] = true;

    while (!bfs.empty()) {
        uint32_t fi = bfs.front();
        bfs.pop();

        const auto& tri = m.faces[fi];
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];

            // Neighbour face should have the reverse directed edge b→a
            auto it = half_edge_face.find(makeEdgeKey(b, a));
            if (it == half_edge_face.end()) {
                // Boundary edge or already flipped; also check same direction
                it = half_edge_face.find(makeEdgeKey(a, b));
                if (it == half_edge_face.end()) continue;

                uint32_t ni = it->second;
                if (ni == fi) continue;
                if (!visited[ni]) {
                    visited[ni] = true;
                    // Same directed edge found → neighbour has wrong winding.
                    // Remove the three OLD forward half-edges before flipping,
                    // then insert the three NEW forward half-edges.
                    auto& ntri = m.faces[ni];
                    for (int ne = 0; ne < 3; ++ne) {
                        half_edge_face.erase(
                            makeEdgeKey(ntri.vi[ne], ntri.vi[(ne + 1) % 3]));
                    }
                    std::swap(ntri.vi[1], ntri.vi[2]);
                    for (int ne = 0; ne < 3; ++ne) {
                        half_edge_face[makeEdgeKey(ntri.vi[ne],
                                                   ntri.vi[(ne + 1) % 3])] = ni;
                    }
                    flipped[ni] = true;
                    bfs.push(ni);
                }
            } else {
                uint32_t ni = it->second;
                if (ni == fi) continue;
                if (!visited[ni]) {
                    visited[ni] = true;
                    bfs.push(ni);
                }
            }
        }
    }

    // Recompute normals for all faces
    for (auto& tri : m.faces) {
        const Vec3d& A = m.vertices[tri.vi[0]].pos;
        const Vec3d& B = m.vertices[tri.vi[1]].pos;
        const Vec3d& C = m.vertices[tri.vi[2]].pos;

        Vec3d ab = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
        Vec3d ac = { C[0]-A[0], C[1]-A[1], C[2]-A[2] };

        Vec3d n = {
            ab[1]*ac[2] - ab[2]*ac[1],
            ab[2]*ac[0] - ab[0]*ac[2],
            ab[0]*ac[1] - ab[1]*ac[0]
        };
        double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (len > 1e-15) {
            tri.normal = { n[0]/len, n[1]/len, n[2]/len };
        }
    }
}

// ---------------------------------------------------------------------------
// fillHoles — fan-triangulate small boundary loops (≤ 20 edges)
// ---------------------------------------------------------------------------
void GeometryEngine::fillHoles(SurfaceManifold& m) {
    if (m.faces.empty()) return;

    // Count edge occurrences (undirected)
    using UEdgeKey = uint64_t;
    auto makeUEdgeKey = [](uint32_t a, uint32_t b) -> UEdgeKey {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    std::unordered_map<UEdgeKey, uint32_t> edge_count;
    edge_count.reserve(m.faces.size() * 6);

    for (const auto& tri : m.faces) {
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            ++edge_count[makeUEdgeKey(a, b)];
        }
    }

    // Collect boundary edges (count == 1) and build adjacency for loop tracing
    // boundary_adj[v] = list of vertices reachable via a boundary edge from v
    std::unordered_map<uint32_t, std::vector<uint32_t>> boundary_adj;
    for (const auto& [key, cnt] : edge_count) {
        if (cnt == 1) {
            uint32_t a = static_cast<uint32_t>(key >> 32);
            uint32_t b = static_cast<uint32_t>(key & 0xFFFFFFFFull);
            boundary_adj[a].push_back(b);
            boundary_adj[b].push_back(a);
        }
    }

    if (boundary_adj.empty()) return;  // already closed

    std::unordered_set<uint32_t> visited_verts;

    // Trace closed boundary loops
    for (auto& [start_v, neighbours] : boundary_adj) {
        if (visited_verts.count(start_v)) continue;

        // Walk the loop
        std::vector<uint32_t> loop;
        loop.push_back(start_v);
        visited_verts.insert(start_v);

        uint32_t prev = start_v;
        uint32_t cur  = neighbours[0];  // pick first neighbour to start

        while (cur != start_v && !visited_verts.count(cur)) {
            loop.push_back(cur);
            visited_verts.insert(cur);

            const auto& adj = boundary_adj[cur];
            uint32_t next = (adj[0] != prev) ? adj[0]
                          : (adj.size() > 1)  ? adj[1]
                          : start_v;
            prev = cur;
            cur  = next;
        }

        // Only fill simple, small holes
        if (loop.size() < 3 || loop.size() > 20) continue;

        // Compute centroid of loop vertices
        Vec3d centroid = {0.0, 0.0, 0.0};
        for (uint32_t vi : loop) {
            const auto& p = m.vertices[vi].pos;
            centroid[0] += p[0];
            centroid[1] += p[1];
            centroid[2] += p[2];
        }
        double inv_n = 1.0 / static_cast<double>(loop.size());
        centroid[0] *= inv_n;
        centroid[1] *= inv_n;
        centroid[2] *= inv_n;

        // Add centroid vertex
        uint32_t cv = static_cast<uint32_t>(m.vertices.size());
        m.vertices.push_back(Vertex{centroid});

        // Fan-triangulate
        uint32_t n = static_cast<uint32_t>(loop.size());
        for (uint32_t i = 0; i < n; ++i) {
            uint32_t a = loop[i];
            uint32_t b = loop[(i + 1) % n];

            Triangle patch_tri;
            patch_tri.vi = { a, b, cv };

            // Compute normal
            const Vec3d& A = m.vertices[a].pos;
            const Vec3d& B = m.vertices[b].pos;
            const Vec3d& C = centroid;

            Vec3d ab = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
            Vec3d ac = { C[0]-A[0], C[1]-A[1], C[2]-A[2] };
            Vec3d nn = {
                ab[1]*ac[2] - ab[2]*ac[1],
                ab[2]*ac[0] - ab[0]*ac[2],
                ab[0]*ac[1] - ab[1]*ac[0]
            };
            double len = std::sqrt(nn[0]*nn[0] + nn[1]*nn[1] + nn[2]*nn[2]);
            if (len > 1e-15) {
                patch_tri.normal = { nn[0]/len, nn[1]/len, nn[2]/len };
            }

            m.faces.push_back(patch_tri);
        }
    }
}

// ---------------------------------------------------------------------------
// checkWatertight — every edge must appear exactly twice
// ---------------------------------------------------------------------------
bool GeometryEngine::checkWatertight(const SurfaceManifold& m) const {
    std::unordered_map<uint64_t, uint32_t> edge_count;
    edge_count.reserve(m.faces.size() * 6);

    auto key = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    for (const auto& tri : m.faces) {
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            ++edge_count[key(a, b)];
        }
    }

    return std::all_of(edge_count.begin(), edge_count.end(),
                       [](const auto& kv) { return kv.second == 2; });
}

// ---------------------------------------------------------------------------
// loadSTL — binary and ASCII
// ---------------------------------------------------------------------------
SurfaceManifold GeometryEngine::loadSTL(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw GeometryError("Cannot open file: " + path);
    }

    // Read first 5 bytes to distinguish ASCII from binary STL
    char header_check[6] = {};
    file.read(header_check, 5);
    if (!file) {
        throw GeometryError("File too short to be a valid STL: " + path);
    }
    file.seekg(0, std::ios::beg);

    bool is_ascii = (std::strncmp(header_check, "solid", 5) == 0);

    // Some binary STLs also start with "solid"; verify by checking file size
    if (is_ascii) {
        file.seekg(0, std::ios::end);
        std::streamsize file_size = file.tellg();
        file.seekg(80, std::ios::beg);

        if (file_size >= 84) {
            uint32_t tri_count = 0;
            file.read(reinterpret_cast<char*>(&tri_count), 4);
            std::streamsize expected = 84 + static_cast<std::streamsize>(tri_count) * 50;
            is_ascii = (file_size != expected) || tri_count == 0;
        }
        file.seekg(0, std::ios::beg);
    }

    SurfaceManifold manifold;

    if (is_ascii) {
        // ASCII STL parser
        std::string line;
        Triangle cur_tri;
        int vertex_idx = 0;
        bool in_facet = false;

        while (std::getline(file, line)) {
            // Trim leading whitespace
            auto pos = line.find_first_not_of(" \t\r\n");
            if (pos == std::string::npos) continue;
            std::string token = line.substr(pos);

            if (token.rfind("facet normal", 0) == 0) {
                in_facet    = true;
                vertex_idx  = 0;
                double nx, ny, nz;
                if (std::sscanf(token.c_str(), "facet normal %lf %lf %lf",
                                &nx, &ny, &nz) == 3) {
                    cur_tri.normal = { nx, ny, nz };
                }
            } else if (in_facet && token.rfind("vertex", 0) == 0) {
                double x, y, z;
                if (std::sscanf(token.c_str(), "vertex %lf %lf %lf",
                                &x, &y, &z) == 3 && vertex_idx < 3) {
                    uint32_t vi = static_cast<uint32_t>(manifold.vertices.size());
                    manifold.vertices.push_back(Vertex{{ x, y, z }});
                    cur_tri.vi[vertex_idx++] = vi;
                }
            } else if (token.rfind("endfacet", 0) == 0 && in_facet) {
                if (vertex_idx == 3) {
                    manifold.faces.push_back(cur_tri);
                }
                in_facet = false;
            }
        }
    } else {
        // Binary STL parser
        // 80-byte header (skip)
        char header[80];
        file.read(header, 80);

        uint32_t tri_count = 0;
        file.read(reinterpret_cast<char*>(&tri_count), 4);

        manifold.vertices.reserve(tri_count * 3);
        manifold.faces.reserve(tri_count);

        for (uint32_t i = 0; i < tri_count; ++i) {
            // 3 floats: normal
            float nx, ny, nz;
            file.read(reinterpret_cast<char*>(&nx), 4);
            file.read(reinterpret_cast<char*>(&ny), 4);
            file.read(reinterpret_cast<char*>(&nz), 4);

            Triangle tri;
            tri.normal = { static_cast<double>(nx),
                           static_cast<double>(ny),
                           static_cast<double>(nz) };

            // 3 × 3 floats: vertices
            for (int v = 0; v < 3; ++v) {
                float vx, vy, vz;
                file.read(reinterpret_cast<char*>(&vx), 4);
                file.read(reinterpret_cast<char*>(&vy), 4);
                file.read(reinterpret_cast<char*>(&vz), 4);

                uint32_t vi = static_cast<uint32_t>(manifold.vertices.size());
                manifold.vertices.push_back(Vertex{{
                    static_cast<double>(vx),
                    static_cast<double>(vy),
                    static_cast<double>(vz)
                }});
                tri.vi[v] = vi;
            }

            // 2-byte attribute byte count (skip)
            uint16_t attr = 0;
            file.read(reinterpret_cast<char*>(&attr), 2);

            if (!file) {
                throw GeometryError("Unexpected EOF reading binary STL at triangle " +
                                    std::to_string(i));
            }

            manifold.faces.push_back(tri);
        }
    }

    if (manifold.faces.empty()) {
        throw GeometryError("STL file contains no triangles: " + path);
    }

    return manifold;
}

// ---------------------------------------------------------------------------
// loadViaGmsh — stub for STEP/IGES (production: link against Gmsh API)
// ---------------------------------------------------------------------------
SurfaceManifold GeometryEngine::loadViaGmsh(const std::string& path) {
    // In production this would call:
    //   gmsh::initialize();
    //   gmsh::model::occ::importShapes(path);
    //   gmsh::model::mesh::generate(2);
    //   ... extract nodes and elements ...
    //   gmsh::finalize();
    //
    // The stub raises a descriptive error so integration is obvious.
    throw GeometryError(
        "STEP/IGES import requires Gmsh or OpenCASCADE linkage. "
        "File: " + path + " — compile with -DUSE_GMSH and link -lgmsh.");
}

// ---------------------------------------------------------------------------
// importCAD
// ---------------------------------------------------------------------------
SurfaceManifold GeometryEngine::importCAD(const std::string& path) {
    emit(progress_, "import", 0.0);

    // Determine format by extension (case-insensitive)
    auto ext_start = path.rfind('.');
    std::string ext;
    if (ext_start != std::string::npos) {
        ext = path.substr(ext_start + 1);
        std::transform(ext.begin(), ext.end(), ext.begin(),
                       [](unsigned char c) { return std::tolower(c); });
    }

    SurfaceManifold manifold;

    if (ext == "stl") {
        emit(progress_, "import:stl", 0.1);
        manifold = loadSTL(path);
    } else if (ext == "step" || ext == "stp" ||
               ext == "iges" || ext == "igs") {
        emit(progress_, "import:cad", 0.1);
        manifold = loadViaGmsh(path);
    } else {
        // Try STL as fallback
        try {
            emit(progress_, "import:stl_fallback", 0.1);
            manifold = loadSTL(path);
        } catch (const GeometryError&) {
            throw GeometryError("Unsupported or unrecognised file format: " + path);
        }
    }

    emit(progress_, "import:done", 1.0);
    return manifold;
}

// ---------------------------------------------------------------------------
// identifyBoundaries
// ---------------------------------------------------------------------------
void GeometryEngine::identifyBoundaries(
        SurfaceManifold& manifold,
        const std::vector<std::pair<std::string, BoundaryType>>& rules) {

    // Helper: case-insensitive substring check
    auto contains_ci = [](const std::string& haystack,
                           const std::string& needle) -> bool {
        if (needle.empty()) return true;
        std::string h = haystack, n = needle;
        std::transform(h.begin(), h.end(), h.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        std::transform(n.begin(), n.end(), n.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return h.find(n) != std::string::npos;
    };

    // Determine BoundaryType from a patch name using keyword heuristics
    auto typeFromName = [&](const std::string& name) -> BoundaryType {
        // Explicit rule table takes priority (first match wins)
        for (const auto& [pattern, btype] : rules) {
            if (contains_ci(name, pattern)) return btype;
        }
        // Built-in keyword heuristics
        if (contains_ci(name, "inlet"))   return BoundaryType::Inlet;
        if (contains_ci(name, "outlet"))  return BoundaryType::Outlet;
        if (contains_ci(name, "sym"))     return BoundaryType::Symmetry;
        if (contains_ci(name, "far"))     return BoundaryType::FarField;
        if (contains_ci(name, "wall"))    return BoundaryType::Wall;
        return BoundaryType::Wall;  // default unlabelled → Wall
    };

    if (manifold.patches.empty()) {
        // No pre-existing patches: assign all faces to a single "wall" patch
        BoundaryPatch wall_patch;
        wall_patch.name = "wall";
        wall_patch.type = BoundaryType::Wall;
        wall_patch.face_ids.resize(manifold.faces.size());
        std::iota(wall_patch.face_ids.begin(), wall_patch.face_ids.end(), 0u);
        manifold.patches.push_back(std::move(wall_patch));
    } else {
        // Update types on existing named patches via rules
        for (auto& patch : manifold.patches) {
            patch.type = typeFromName(patch.name);
        }
    }
}

// ---------------------------------------------------------------------------
// repairManifold
// ---------------------------------------------------------------------------
bool GeometryEngine::repairManifold(SurfaceManifold& manifold) {
    emit(progress_, "repair:merge_vertices", 0.0);
    mergeVertices(manifold);

    emit(progress_, "repair:fix_winding", 0.33);
    fixWinding(manifold);

    emit(progress_, "repair:fill_holes", 0.66);
    fillHoles(manifold);

    emit(progress_, "repair:check", 0.9);
    manifold.is_watertight = checkWatertight(manifold);

    emit(progress_, "repair:done", 1.0);
    auto report = analyzeManifold(manifold);
    if (report.non_manifold_edges > 0) {
        emit(progress_, "repair:warning_non_manifold", 0.95);
    }
    return manifold.is_watertight;
}

// ---------------------------------------------------------------------------
// patchCentroid — area-weighted centroid
// ---------------------------------------------------------------------------
Vec3d GeometryEngine::patchCentroid(const SurfaceManifold& m,
                                    const BoundaryPatch& p) const {
    if (p.face_ids.empty()) return {0.0, 0.0, 0.0};

    Vec3d  weighted_sum = {0.0, 0.0, 0.0};
    double total_area   = 0.0;

    for (uint32_t fi : p.face_ids) {
        if (fi >= static_cast<uint32_t>(m.faces.size())) continue;

        const auto& tri = m.faces[fi];
        const Vec3d& A  = m.vertices[tri.vi[0]].pos;
        const Vec3d& B  = m.vertices[tri.vi[1]].pos;
        const Vec3d& C  = m.vertices[tri.vi[2]].pos;

        // Triangle centroid
        Vec3d fc = {
            (A[0] + B[0] + C[0]) / 3.0,
            (A[1] + B[1] + C[1]) / 3.0,
            (A[2] + B[2] + C[2]) / 3.0
        };

        // Triangle area = 0.5 * |AB × AC|
        Vec3d ab = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
        Vec3d ac = { C[0]-A[0], C[1]-A[1], C[2]-A[2] };
        Vec3d cross = {
            ab[1]*ac[2] - ab[2]*ac[1],
            ab[2]*ac[0] - ab[0]*ac[2],
            ab[0]*ac[1] - ab[1]*ac[0]
        };
        double area = 0.5 * std::sqrt(cross[0]*cross[0] +
                                      cross[1]*cross[1] +
                                      cross[2]*cross[2]);

        weighted_sum[0] += area * fc[0];
        weighted_sum[1] += area * fc[1];
        weighted_sum[2] += area * fc[2];
        total_area += area;
    }

    if (total_area < 1e-15) return {0.0, 0.0, 0.0};

    return {
        weighted_sum[0] / total_area,
        weighted_sum[1] / total_area,
        weighted_sum[2] / total_area
    };
}

// ---------------------------------------------------------------------------
// simplifyLoop — remove collinear midpoints from an ordered vertex loop
// ---------------------------------------------------------------------------
static std::vector<uint32_t> simplifyLoop(
        const std::vector<uint32_t>& loop,
        const SurfaceManifold& m,
        double eps = 1e-6) {
    if (loop.size() <= 3) return loop;

    const size_t n = loop.size();
    std::vector<uint32_t> result;
    result.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        const Vec3d& A = m.vertices[loop[(i + n - 1) % n]].pos;
        const Vec3d& B = m.vertices[loop[i]].pos;
        const Vec3d& C = m.vertices[loop[(i + 1) % n]].pos;

        Vec3d ab = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
        Vec3d ac = { C[0]-A[0], C[1]-A[1], C[2]-A[2] };
        double cx = ab[1]*ac[2] - ab[2]*ac[1];
        double cy = ab[2]*ac[0] - ab[0]*ac[2];
        double cz = ab[0]*ac[1] - ab[1]*ac[0];
        double cross2 = cx*cx + cy*cy + cz*cz;
        double ab2    = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
        double ac2    = ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2];

        // Keep vertex if it is NOT collinear (|cross| >= eps * |AB| * |AC|)
        if (cross2 >= eps * eps * ab2 * ac2) {
            result.push_back(loop[i]);
        }
    }

    return (result.size() >= 3) ? result : loop;
}

// ---------------------------------------------------------------------------
// identifyFlatFaces — BFS flat-face grouping via edge-adjacency
// ---------------------------------------------------------------------------
std::vector<std::vector<uint32_t>> GeometryEngine::identifyFlatFaces(
        const SurfaceManifold& m, double epsilon) const {

    const uint32_t nf = static_cast<uint32_t>(m.faces.size());
    if (nf == 0) return {};

    // Undirected edge (sorted vertex pair) → list of face indices
    auto edgeKey = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    std::unordered_map<uint64_t, std::vector<uint32_t>> edgeToFaces;
    edgeToFaces.reserve(nf * 4);

    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& tri = m.faces[fi];
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            edgeToFaces[edgeKey(a, b)].push_back(fi);
        }
    }

    const double cos_thresh = 1.0 - epsilon;

    auto dot3 = [](const Vec3d& a, const Vec3d& b) -> double {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };

    std::vector<bool> visited(nf, false);
    std::vector<std::vector<uint32_t>> groups;

    for (uint32_t i = 0; i < nf; ++i) {
        if (visited[i]) continue;

        std::vector<uint32_t> group;
        std::vector<uint32_t> queue;
        queue.push_back(i);
        visited[i] = true;

        uint32_t head = 0;
        while (head < static_cast<uint32_t>(queue.size())) {
            uint32_t fi = queue[head++];
            group.push_back(fi);

            const auto& tri = m.faces[fi];
            for (int e = 0; e < 3; ++e) {
                uint32_t a = tri.vi[e];
                uint32_t b = tri.vi[(e + 1) % 3];
                auto it = edgeToFaces.find(edgeKey(a, b));
                if (it == edgeToFaces.end()) continue;

                for (uint32_t ni : it->second) {
                    if (visited[ni]) continue;
                    if (dot3(tri.normal, m.faces[ni].normal) >= cos_thresh) {
                        visited[ni] = true;
                        queue.push_back(ni);
                    }
                }
            }
        }

        groups.push_back(std::move(group));
    }

    return groups;
}

// ---------------------------------------------------------------------------
// extractAllBoundaryLoops — one-pass boundary loop extraction for all groups
// ---------------------------------------------------------------------------
std::vector<std::vector<std::vector<uint32_t>>>
GeometryEngine::extractAllBoundaryLoops(
        const SurfaceManifold& m,
        const std::vector<std::vector<uint32_t>>& groups,
        bool simplify) const {

    const uint32_t nf = static_cast<uint32_t>(m.faces.size());
    const uint32_t ng = static_cast<uint32_t>(groups.size());

    if (nf == 0 || ng == 0) return {};

    // Map each triangle index → its group index
    std::vector<uint32_t> triToGroup(nf, UINT32_MAX);
    for (uint32_t gi = 0; gi < ng; ++gi) {
        for (uint32_t fi : groups[gi]) {
            if (fi < nf) triToGroup[fi] = gi;
        }
    }

    // Build directed-edge → owning face index (a→b packed into uint64_t)
    auto dKey = [](uint32_t a, uint32_t b) -> uint64_t {
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    std::unordered_map<uint64_t, uint32_t> dirEdgeFace;
    dirEdgeFace.reserve(nf * 6);
    for (uint32_t fi = 0; fi < nf; ++fi) {
        const auto& tri = m.faces[fi];
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            dirEdgeFace[dKey(a, b)] = fi;
        }
    }

    // For each group, build a directed adjacency map of boundary edges (v → next v).
    // A directed edge a→b is a boundary edge of group gi when its reverse b→a
    // either does not exist or belongs to a different group.
    std::vector<std::unordered_map<uint32_t, uint32_t>> groupNext(ng);

    for (const auto& kv : dirEdgeFace) {
        uint32_t a  = static_cast<uint32_t>(kv.first >> 32);
        uint32_t b  = static_cast<uint32_t>(kv.first & 0xFFFFFFFFull);
        uint32_t fi = kv.second;
        uint32_t gi = triToGroup[fi];
        if (gi == UINT32_MAX) continue;

        auto rev = dirEdgeFace.find(dKey(b, a));
        bool interior = (rev != dirEdgeFace.end() && triToGroup[rev->second] == gi);
        if (!interior) {
            groupNext[gi][a] = b; // boundary: a → b for this group
        }
    }

    // Chain boundary edges into ordered loops
    std::vector<std::vector<std::vector<uint32_t>>> allLoops(ng);

    for (uint32_t gi = 0; gi < ng; ++gi) {
        auto& adj = groupNext[gi];
        std::unordered_set<uint32_t> visited;

        for (auto& startKV : adj) {
            uint32_t start = startKV.first;
            if (visited.count(start)) continue;

            std::vector<uint32_t> loop;
            uint32_t curr = start;

            while (!visited.count(curr)) {
                auto it = adj.find(curr);
                if (it == adj.end()) break;
                visited.insert(curr);
                loop.push_back(curr);
                curr = it->second;
            }

            if (loop.size() >= 3) {
                if (simplify) loop = simplifyLoop(loop, m);
                allLoops[gi].push_back(std::move(loop));
            }
        }
    }

    return allLoops;
}

// ---------------------------------------------------------------------------
// analyzeManifold — classify every undirected edge by valence
// ---------------------------------------------------------------------------
GeometryEngine::ManifoldReport GeometryEngine::analyzeManifold(
        const SurfaceManifold& m) const {

    // Build undirected edge → occurrence count
    auto makeUEdgeKey = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    std::unordered_map<uint64_t, int> edge_count;
    edge_count.reserve(m.faces.size() * 6);

    for (const auto& tri : m.faces) {
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            ++edge_count[makeUEdgeKey(a, b)];
        }
    }

    int boundary_edges    = 0;
    int non_manifold_edges = 0;
    for (const auto& [key, cnt] : edge_count) {
        if (cnt == 1)       ++boundary_edges;
        else if (cnt >= 3)  ++non_manifold_edges;
    }

    // Count isolated vertices: not referenced by any face
    const uint32_t nv = static_cast<uint32_t>(m.vertices.size());
    std::vector<bool> used(nv, false);
    for (const auto& tri : m.faces) {
        for (auto vi : tri.vi) {
            if (vi < nv) used[vi] = true;
        }
    }
    int isolated_vertices = static_cast<int>(
        std::count(used.begin(), used.end(), false));

    return ManifoldReport{boundary_edges, non_manifold_edges, isolated_vertices};
}

// ---------------------------------------------------------------------------
// nonManifoldEdges — return every edge whose valence != 2
// ---------------------------------------------------------------------------
std::vector<Edge> GeometryEngine::nonManifoldEdges(
        const SurfaceManifold& m) const {

    auto makeUEdgeKey = [](uint32_t a, uint32_t b) -> uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
    };

    std::unordered_map<uint64_t, int> edge_count;
    edge_count.reserve(m.faces.size() * 6);

    for (const auto& tri : m.faces) {
        for (int e = 0; e < 3; ++e) {
            uint32_t a = tri.vi[e];
            uint32_t b = tri.vi[(e + 1) % 3];
            ++edge_count[makeUEdgeKey(a, b)];
        }
    }

    std::vector<Edge> result;
    for (const auto& [key, cnt] : edge_count) {
        if (cnt != 2) {
            uint32_t a = static_cast<uint32_t>(key >> 32);
            uint32_t b = static_cast<uint32_t>(key & 0xFFFFFFFFull);
            result.push_back(Edge{{a, b}});
        }
    }
    return result;
}

} // namespace cfd
