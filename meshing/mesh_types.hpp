#pragma once
#include "../geometry/types.hpp"
#include <vector>
#include <cstdint>

namespace cfd {

// ---------------------------------------------------------------------------
// PrismElement — a 6-node wedge element produced by boundary layer inflation.
//
//   Bottom face (wall side):  nodes[0], nodes[1], nodes[2]
//   Top face    (flow side):  nodes[3], nodes[4], nodes[5]
//
//   The triangle (nodes[i], nodes[i+3]) is a lateral quad face of the prism.
// ---------------------------------------------------------------------------
struct PrismElement {
    std::array<uint32_t, 6> nodes;
    uint32_t layer;        // 0 = first layer (nearest wall)
    uint32_t wall_face_id; // index of the source triangle on the wall surface
};

// ---------------------------------------------------------------------------
// TetElement — a 4-node tetrahedron used in the core volume mesh.
// ---------------------------------------------------------------------------
struct TetElement {
    std::array<uint32_t, 4> nodes;
};

// ---------------------------------------------------------------------------
// PrismMesh — all boundary-layer data produced by inflatePrisms().
//
//   vertices        : every node (original wall nodes first, then inflated).
//   wall_node_ids   : indices (into vertices) of the original wall surface nodes.
//   outer_node_ids  : indices of the outermost inflated layer — these become
//                     the inner boundary for the tet core.
// ---------------------------------------------------------------------------
struct PrismMesh {
    std::vector<Vertex>       vertices;
    std::vector<PrismElement> prisms;
    std::vector<uint32_t>     wall_node_ids;   // original wall surface nodes
    std::vector<uint32_t>     outer_node_ids;  // outermost inflated layer nodes
    int    num_layers;
    double first_cell_height;
    double growth_ratio;
};

// ---------------------------------------------------------------------------
// VolumeMesh — the complete hybrid (prism + tet) volume mesh.
// ---------------------------------------------------------------------------
struct VolumeMesh {
    std::vector<Vertex>       vertices;
    std::vector<TetElement>   tets;
    std::vector<PrismElement> prisms;   // boundary layer prisms (copied from PrismMesh)
    std::vector<BoundaryPatch> patches; // inherited from surface + updated node ids
    // Outer domain bounding box used by generateTetCore (set after construction).
    // Faces whose vertices all lie on any face of this box are domain walls.
    double domain_lo[3] = {0,0,0};
    double domain_hi[3] = {0,0,0};
};

// ---------------------------------------------------------------------------
// BLParameters — derived boundary layer sizing parameters.
//
//   total_thickness = h1 * (r^n - 1) / (r - 1)
// ---------------------------------------------------------------------------
struct BLParameters {
    double first_cell_height; // h1 [m]
    double growth_ratio;      // r, e.g. 1.2
    int    num_layers;
    double total_thickness;   // cumulative BL thickness [m]
};

} // namespace cfd
