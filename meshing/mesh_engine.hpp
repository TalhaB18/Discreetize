#pragma once
#include "mesh_types.hpp"
#include <functional>
#include <stdexcept>
#include <string>

namespace cfd {

// ---------------------------------------------------------------------------
// MeshEngine — Boundary Layer inflation + Tet core generation.
//
// Typical call sequence:
//   MeshEngine eng([](auto msg, auto pct){ /* log */ });
//   auto params = MeshEngine::computeBLParameters(cfg);
//   auto bl     = eng.inflatePrisms(surface, params);
//   auto vol    = eng.generateTetCore(bl, cfg.core_mesh_size);
//   auto report = eng.checkQuality(vol);
// ---------------------------------------------------------------------------
class MeshEngine {
public:
    // Progress callback: (stage_description, fraction_complete [0..1])
    using ProgressCb = std::function<void(const std::string&, double)>;

    explicit MeshEngine(ProgressCb cb = nullptr);

    // -----------------------------------------------------------------------
    // y+ physics: first cell height from flow conditions.
    //
    //   C_f    = 0.0592 * Re^(-0.2)              (Schlichting turbulent flat plate)
    //   tau_w  = 0.5 * rho * U_inf^2 * C_f        rho = 1.225 kg/m^3
    //   u_tau  = sqrt(tau_w / rho)                (friction velocity)
    //   h1     = y+ * nu / u_tau
    //
    // Parameters:
    //   Re    — Reynolds number  U_inf * L / nu  (dimensionless)
    //   yplus — target y+ (< 1 for viscous sublayer, 30–300 for wall functions)
    //   nu    — kinematic viscosity [m^2/s]
    //   U_inf — freestream velocity [m/s]
    //
    // Returns h1 [m].
    // -----------------------------------------------------------------------
    [[nodiscard]] static double computeFirstCellHeight(double Re, double yplus,
                                                       double nu, double U_inf);

    // Derive the full BLParameters set from a PipelineConfig.
    [[nodiscard]] static BLParameters computeBLParameters(const PipelineConfig& cfg);

    // -----------------------------------------------------------------------
    // Prism inflation — extrude every wall triangle into N prism layers.
    //
    // For each wall face triangle (v0, v1, v2) with area-weighted vertex normal n:
    //   Cumulative offset at layer k:  h_k = h1 * (r^k - 1) / (r - 1)
    //   Node position at layer k:      pos_k = pos_wall + h_k * n
    //
    // Nodes shared between adjacent wall faces are deduplicated via a
    // position → index map (tolerance 1e-12 m).
    //
    // After extrusion, smoothPrismLayers() is called to repair tangling.
    // -----------------------------------------------------------------------
    PrismMesh inflatePrisms(const SurfaceManifold& surface,
                            const BLParameters& params);

    // -----------------------------------------------------------------------
    // Tetrahedral core — fill the interior with tets.
    //
    // Seed points are placed on a regular Cartesian grid with spacing
    // target_size, clipped to the bounding box of the outer prism shell
    // and filtered to exclude points already inside the prism layer.
    // Bowyer–Watson incremental Delaunay connects all outer-shell nodes
    // plus interior seeds into a conforming tet mesh.
    // -----------------------------------------------------------------------
    VolumeMesh generateTetCore(const PrismMesh& bl_mesh, double target_size);

    // -----------------------------------------------------------------------
    // Quality assessment
    // -----------------------------------------------------------------------
    struct QualityReport {
        double max_skewness;       // 0 = equilateral, 1 = degenerate
        double avg_skewness;
        double max_aspect_ratio;   // longest edge / shortest altitude
        double avg_aspect_ratio;
        double min_orthogonality;  // degrees; > 70° is acceptable
        int    bad_cells;          // count of cells with skewness > 0.85
    };

    [[nodiscard]] QualityReport checkQuality(const VolumeMesh& mesh) const;

private:
    ProgressCb progress_;

    // Emit a progress notification (no-op when progress_ is null).
    void notify(const std::string& stage, double fraction) const;

    // Area-weighted vertex normals, considering only the supplied face indices.
    std::vector<Vec3d> computeVertexNormals(
        const SurfaceManifold& surface,
        const std::vector<uint32_t>& wall_face_ids) const;

    // Laplacian smoothing of inflated positions to reduce layer tangling.
    // Only moves non-wall nodes; wall nodes are pinned.
    void smoothPrismLayers(PrismMesh& mesh, int iterations = 3) const;

    // Bowyer–Watson incremental insertion into an existing tet set.
    // pts must be appended to vol.vertices BEFORE this call.
    void bowyerWatson(VolumeMesh& vol, const std::vector<Vec3d>& interior_pts) const;

    // Per-element metrics (static — no mesh context needed)
    static double tetSkewness(const std::array<Vec3d, 4>& verts);
    static double tetAspectRatio(const std::array<Vec3d, 4>& verts);

    // Vector arithmetic helpers
    static Vec3d cross(const Vec3d& a, const Vec3d& b);
    static double dot(const Vec3d& a, const Vec3d& b);
    static double norm(const Vec3d& a);
    static Vec3d normalize(const Vec3d& a);
    static Vec3d add(const Vec3d& a, const Vec3d& b);
    static Vec3d sub(const Vec3d& a, const Vec3d& b);
    static Vec3d scale(const Vec3d& a, double s);
    static Vec3d midpoint(const Vec3d& a, const Vec3d& b);

public:
    // Circumsphere of a tet: returns {centre, radius}
    static std::pair<Vec3d, double> circumsphere(const std::array<Vec3d, 4>& v);
};

// ---------------------------------------------------------------------------
// MeshError — thrown for unrecoverable meshing failures.
// ---------------------------------------------------------------------------
class MeshError : public std::runtime_error {
public:
    explicit MeshError(const std::string& m) : std::runtime_error(m) {}
};

} // namespace cfd
