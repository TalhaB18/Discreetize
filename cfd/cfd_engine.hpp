#pragma once
#include "../meshing/mesh_types.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <variant>
#include <functional>

namespace cfd {

// Dirichlet: fixed value.  Neumann: fixed gradient (usually 0 = zero-flux)
struct DirichletBC { double value; };
struct NeumannBC   { double gradient; };  // dφ/dn = gradient
using BoundaryCondition = std::variant<DirichletBC, NeumannBC>;

struct FieldBC {
    std::string field_name;   // "U", "p", "k", "epsilon", "omega"
    std::unordered_map<std::string, BoundaryCondition> patch_bcs;
    // patch_name -> BC
};

// FVM cell: stores centroid, volume, face connectivity
struct FVMCell {
    Vec3d    centroid;
    double   volume;
    std::vector<uint32_t> face_ids;
    bool     is_prism;   // false = tet
};

struct FVMFace {
    Vec3d    centroid;
    Vec3d    normal;     // outward from owner cell
    double   area;
    uint32_t owner;      // cell index
    int32_t  neighbour;  // -1 if boundary face
    std::string patch_name; // set if boundary face
};

struct CFDMesh {
    std::vector<Vertex>   points;
    std::vector<FVMCell>  cells;
    std::vector<FVMFace>  faces;
    std::vector<BoundaryPatch> boundary_patches;
    std::vector<FieldBC>  field_bcs;
};

class CFDEngine {
public:
    using ProgressCb = std::function<void(const std::string&, double)>;
    explicit CFDEngine(ProgressCb cb = nullptr);

    // Convert VolumeMesh (tets+prisms) into FVM data structure
    CFDMesh mapMesh(const VolumeMesh& vol);

    // Apply standard boundary conditions for external aerodynamics:
    //   Wall:     U=Dirichlet(0), p=Neumann(0)
    //   Inlet:    U=Dirichlet(U_inf), p=Neumann(0)
    //   Outlet:   U=Neumann(0), p=Dirichlet(0)
    //   Symmetry: all fields Neumann(0)
    void applyBoundaryConditions(CFDMesh& mesh,
                                 double U_inf,
                                 const std::string& turbulence_model = "kOmegaSST");

    // Export in requested format
    void exportMesh(const CFDMesh& mesh,
                    const std::string& path,
                    const std::string& format);  // "openfoam", "su2", "fluent"

    // Validate: check face normals consistent, no negative volumes
    bool validateMesh(const CFDMesh& mesh) const;

private:
    ProgressCb progress_;

    void exportOpenFOAM(const CFDMesh& mesh, const std::string& dir) const;
    void exportSU2(const CFDMesh& mesh, const std::string& path) const;
    void exportFluent(const CFDMesh& mesh, const std::string& path) const;

    // Compute tet/prism centroid and volume
    static Vec3d tetCentroid(const std::array<Vec3d,4>& v);
    static double tetVolume(const std::array<Vec3d,4>& v);
    static Vec3d prismCentroid(const std::array<Vec3d,6>& v);
    static double prismVolume(const std::array<Vec3d,6>& v);
    static Vec3d faceNormal(const Vec3d& a, const Vec3d& b, const Vec3d& c);
};

class CFDError : public std::runtime_error {
public:
    explicit CFDError(const std::string& m) : std::runtime_error(m) {}
};

// -----------------------------------------------------------------------
// FVM Solver hook — abstract interface for plugging in external solvers
// -----------------------------------------------------------------------
class FVMSolverHook {
public:
    virtual ~FVMSolverHook() = default;

    // Called once before time-stepping with the CFD mesh
    virtual void initialize(const CFDMesh& mesh) = 0;

    // Solve one pseudo-time step; returns L2 residual
    virtual double step(int iteration) = 0;

    // Extract field values at cell centres after convergence
    virtual std::vector<double> getField(const std::string& name) const = 0;

    // Convergence check
    virtual bool isConverged(double tol = 1e-6) const = 0;
};

// Concrete stub solver — simple Laplace for pressure (demo only)
class LaplaceSolverHook : public FVMSolverHook {
public:
    void initialize(const CFDMesh& mesh) override;
    double step(int iteration) override;
    std::vector<double> getField(const std::string& name) const override;
    bool isConverged(double tol = 1e-6) const override;
private:
    const CFDMesh* mesh_ = nullptr;
    std::vector<double> pressure_;
    double residual_ = 1.0;
};

} // namespace cfd
