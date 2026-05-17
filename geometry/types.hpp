#pragma once
#include <array>
#include <vector>
#include <string>
#include <cstdint>

namespace cfd {

using Vec3d = std::array<double, 3>;
using Vec3i = std::array<int32_t, 3>;

struct Vertex   { Vec3d pos; };
struct Triangle { std::array<uint32_t, 3> vi; Vec3d normal; };
struct Edge     { std::array<uint32_t, 2> vi; };

enum class BoundaryType {
    Wall,        // no-slip — boundary layer required
    Inlet,       // Dirichlet velocity
    Outlet,      // Neumann pressure
    Symmetry,    // zero-gradient normal
    FarField,    // free-stream
    Unknown
};

struct BoundaryPatch {
    std::string      name;
    BoundaryType     type;
    std::vector<uint32_t> face_ids;  // indices into SurfaceManifold::faces
};

struct SurfaceManifold {
    std::vector<Vertex>   vertices;
    std::vector<Triangle> faces;
    std::vector<BoundaryPatch> patches;
    bool is_watertight = false;
};

struct PipelineConfig {
    std::string input_file;       // path to STEP/IGES/STL
    double      reynolds_number;
    double      target_yplus;     // typically <1 for viscous sublayer
    double      kinematic_viscosity; // nu [m^2/s]
    double      freestream_velocity; // U_inf [m/s]
    int         bl_layers;        // number of prism layers
    double      bl_growth_ratio;  // e.g. 1.2
    double      core_mesh_size;   // target tet size
    std::string output_format;    // "openfoam", "fluent", "su2"
};

} // namespace cfd
