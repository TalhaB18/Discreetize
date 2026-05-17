#pragma once
#include "../geometry/types.hpp"
#include <vector>
#include <array>
#include <string>

namespace cfd {

struct QuadFace { std::array<uint32_t, 4> vi; };

struct QuadMesh {
    std::vector<Vertex>   vertices;
    std::vector<QuadFace> quads;
};

// Convert a triangle surface mesh to an all-quad mesh via centroid subdivision.
// Each triangle (v0,v1,v2) becomes 3 quads using centroid c and edge midpoints:
//   (v0, m01, c, m20), (v1, m12, c, m01), (v2, m20, c, m12)
class QuadMesher {
public:
    QuadMesh convert(const SurfaceManifold& m) const;
    void writeVTU(const QuadMesh& qm, const std::string& path) const;
};

} // namespace cfd
