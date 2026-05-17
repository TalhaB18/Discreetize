#pragma once
#include "types.hpp"
#include <memory>
#include <functional>
#include <stdexcept>

namespace cfd {

class GeometryEngine {
public:
    // Progress callback: (stage_name, 0.0-1.0)
    using ProgressCb = std::function<void(const std::string&, double)>;

    explicit GeometryEngine(ProgressCb cb = nullptr);
    ~GeometryEngine();

    // Import from file. Supports STL (native), STEP/IGES (via stub or OCC).
    // Returns watertight SurfaceManifold or throws GeometryError.
    SurfaceManifold importCAD(const std::string& path);

    // Tag surface patches by name patterns or normal direction heuristics.
    // walls: faces with inward normals matching "wall_*" or all unlabelled faces
    void identifyBoundaries(SurfaceManifold& manifold,
                            const std::vector<std::pair<std::string,BoundaryType>>& rules);

    // Attempt to close holes, merge duplicate vertices, fix winding.
    // Returns true if manifold is watertight after repair.
    [[nodiscard]] bool repairManifold(SurfaceManifold& manifold);

    // Compute surface area weighted centroid per patch (useful for BC placement)
    Vec3d patchCentroid(const SurfaceManifold& m, const BoundaryPatch& p) const;

    // Merge duplicate vertices (required before identifyFlatFaces on raw STL data).
    void mergeVertices(SurfaceManifold& m, double tol = 1e-8);

    // BFS flood-fill grouping: adjacent triangles with nearly-parallel normals
    // (dot product >= 1 - epsilon) form one flat face group.
    // Returns one vector of triangle indices per group.
    // Call mergeVertices first so that edges are shared between triangles.
    [[nodiscard]] std::vector<std::vector<uint32_t>> identifyFlatFaces(
        const SurfaceManifold& m, double epsilon = 0.001) const;

    // For each face group, extract ordered boundary vertex-index loops.
    // Returns allLoops[group][loop_index] = ordered list of vertex indices.
    // Inner loops (holes) are included as additional loop entries.
    // Collinear midpoints are removed when simplify=true.
    // Efficient single-pass over edges (O(total_triangles)).
    std::vector<std::vector<std::vector<uint32_t>>> extractAllBoundaryLoops(
        const SurfaceManifold& m,
        const std::vector<std::vector<uint32_t>>& groups,
        bool simplify = true) const;

    struct ManifoldReport {
        int boundary_edges;       // edges shared by exactly 1 triangle
        int non_manifold_edges;   // edges shared by 3+ triangles
        int isolated_vertices;    // vertices referenced by no face
        bool is_manifold() const { return non_manifold_edges == 0 && boundary_edges == 0; }
    };

    ManifoldReport analyzeManifold(const SurfaceManifold& m) const;
    // Returns detailed report without modifying m.

    std::vector<Edge> nonManifoldEdges(const SurfaceManifold& m) const;
    // Returns the list of edges shared by != 2 triangles.
    // Empty for a valid closed manifold.

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    ProgressCb progress_;

    SurfaceManifold loadSTL(const std::string& path);
    SurfaceManifold loadViaGmsh(const std::string& path);
    void fixWinding(SurfaceManifold& m);
    void fillHoles(SurfaceManifold& m);
    bool checkWatertight(const SurfaceManifold& m) const;
};

class GeometryError : public std::runtime_error {
public:
    explicit GeometryError(const std::string& msg) : std::runtime_error(msg) {}
};

} // namespace cfd
