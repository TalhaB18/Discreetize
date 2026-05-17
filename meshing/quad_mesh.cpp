#include "quad_mesh.hpp"
#include <fstream>
#include <map>
#include <cstdio>

namespace cfd {

QuadMesh QuadMesher::convert(const SurfaceManifold& m) const
{
    QuadMesh qm;
    // Pre-allocate: 3 quads per triangle, each using 5 new vertices
    // (1 centroid + 3 midpoints, midpoints shared between triangles)
    // Use a map to deduplicate edge midpoints.
    qm.vertices = m.vertices;  // copy original vertices

    // midpoint map: {min_vi, max_vi} -> midpoint vertex index
    std::map<std::pair<uint32_t,uint32_t>, uint32_t> midMap;

    auto getMid = [&](uint32_t a, uint32_t b) -> uint32_t {
        auto key = (a < b) ? std::make_pair(a,b) : std::make_pair(b,a);
        auto it = midMap.find(key);
        if(it != midMap.end()) return it->second;
        const Vec3d& pa = qm.vertices[a].pos;
        const Vec3d& pb = qm.vertices[b].pos;
        Vertex v;
        v.pos = {(pa[0]+pb[0])*0.5, (pa[1]+pb[1])*0.5, (pa[2]+pb[2])*0.5};
        uint32_t idx = static_cast<uint32_t>(qm.vertices.size());
        qm.vertices.push_back(v);
        midMap[key] = idx;
        return idx;
    };

    qm.quads.reserve(m.faces.size() * 3);

    for(const auto& tri : m.faces){
        uint32_t v0 = tri.vi[0], v1 = tri.vi[1], v2 = tri.vi[2];

        // Centroid
        const Vec3d& p0 = qm.vertices[v0].pos;
        const Vec3d& p1 = qm.vertices[v1].pos;
        const Vec3d& p2 = qm.vertices[v2].pos;
        Vertex vc;
        vc.pos = {(p0[0]+p1[0]+p2[0])/3.0,
                  (p0[1]+p1[1]+p2[1])/3.0,
                  (p0[2]+p1[2]+p2[2])/3.0};
        uint32_t c = static_cast<uint32_t>(qm.vertices.size());
        qm.vertices.push_back(vc);

        uint32_t m01 = getMid(v0, v1);
        uint32_t m12 = getMid(v1, v2);
        uint32_t m20 = getMid(v2, v0);

        qm.quads.push_back({v0, m01, c, m20});
        qm.quads.push_back({v1, m12, c, m01});
        qm.quads.push_back({v2, m20, c, m12});
    }

    return qm;
}

void QuadMesher::writeVTU(const QuadMesh& qm, const std::string& path) const
{
    std::ofstream f(path);
    if(!f) return;

    const size_t nv = qm.vertices.size();
    const size_t nq = qm.quads.size();

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << nv
      << "\" NumberOfCells=\"" << nq << "\">\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for(const auto& v : qm.vertices)
        f << "          " << v.pos[0] << " " << v.pos[1] << " " << v.pos[2] << "\n";

    f << "        </DataArray>\n"
      << "      </Points>\n"
      << "      <Cells>\n"
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";

    for(const auto& q : qm.quads)
        f << "          " << q.vi[0] << " " << q.vi[1] << " " << q.vi[2] << " " << q.vi[3] << "\n";

    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
    for(size_t i=1; i<=nq; i++) f << (i*4) << " ";
    f << "\n        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
    for(size_t i=0; i<nq; i++) f << "9 ";  // VTK_QUAD = 9
    f << "\n        </DataArray>\n"
      << "      </Cells>\n"
      << "    </Piece>\n"
      << "  </UnstructuredGrid>\n"
      << "</VTKFile>\n";
}

} // namespace cfd
