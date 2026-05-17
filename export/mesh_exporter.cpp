#include "mesh_exporter.hpp"
#include "../analysis/curvature.hpp"
#include <fstream>
#include <stdexcept>
#include <cstring>
#include <cstdint>

namespace cfd {

// ---------------------------------------------------------------------------
// STLExporter
// ---------------------------------------------------------------------------

static_assert(sizeof(float) == 4, "float must be 32-bit");

void STLExporter::write(const SurfaceManifold& m, const std::string& path) const {
    if (mode_ == Mode::Binary)
        writeBinary(m, path);
    else
        writeASCII(m, path);
}

std::string STLExporter::extension() const {
    return ".stl";
}

void STLExporter::writeBinary(const SurfaceManifold& m, const std::string& path) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("STLExporter: cannot open " + path);

    char header[80]{};
    std::strncpy(header, "MeshMaster STL", sizeof(header) - 1);
    f.write(header, 80);

    auto count = static_cast<uint32_t>(m.faces.size());
    f.write(reinterpret_cast<const char*>(&count), sizeof(count));

    for (const auto& tri : m.faces) {
        float nx = static_cast<float>(tri.normal[0]);
        float ny = static_cast<float>(tri.normal[1]);
        float nz = static_cast<float>(tri.normal[2]);
        f.write(reinterpret_cast<const char*>(&nx), 4);
        f.write(reinterpret_cast<const char*>(&ny), 4);
        f.write(reinterpret_cast<const char*>(&nz), 4);

        for (int v = 0; v < 3; ++v) {
            const auto& pos = m.vertices[tri.vi[v]].pos;
            float vx = static_cast<float>(pos[0]);
            float vy = static_cast<float>(pos[1]);
            float vz = static_cast<float>(pos[2]);
            f.write(reinterpret_cast<const char*>(&vx), 4);
            f.write(reinterpret_cast<const char*>(&vy), 4);
            f.write(reinterpret_cast<const char*>(&vz), 4);
        }

        uint16_t attr = 0;
        f.write(reinterpret_cast<const char*>(&attr), sizeof(attr));
    }

    if (!f) throw std::runtime_error("STLExporter: write error on " + path);
}

void STLExporter::writeASCII(const SurfaceManifold& m, const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("STLExporter: cannot open " + path);

    f << "solid meshmaster\n";
    for (const auto& tri : m.faces) {
        f << "facet normal "
          << tri.normal[0] << " "
          << tri.normal[1] << " "
          << tri.normal[2] << "\n";
        f << "  outer loop\n";
        for (int v = 0; v < 3; ++v) {
            const auto& pos = m.vertices[tri.vi[v]].pos;
            f << "    vertex " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }
        f << "  endloop\n";
        f << "endfacet\n";
    }
    f << "endsolid meshmaster\n";

    if (!f) throw std::runtime_error("STLExporter: write error on " + path);
}

// ---------------------------------------------------------------------------
// OBJExporter
// ---------------------------------------------------------------------------

void OBJExporter::write(const SurfaceManifold& m, const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("OBJExporter: cannot open " + path);

    for (const auto& vert : m.vertices) {
        f << "v " << vert.pos[0] << " " << vert.pos[1] << " " << vert.pos[2] << "\n";
    }

    for (const auto& tri : m.faces) {
        f << "vn " << tri.normal[0] << " " << tri.normal[1] << " " << tri.normal[2] << "\n";
    }

    if (groups_ && !groups_->empty()) {
        for (std::size_t g = 0; g < groups_->size(); ++g) {
            f << "g group_" << g << "\n";
            for (uint32_t fi : (*groups_)[g]) {
                const auto& tri = m.faces[fi];
                uint32_t ni = fi + 1;
                f << "f "
                  << (tri.vi[0] + 1) << "//" << ni << " "
                  << (tri.vi[1] + 1) << "//" << ni << " "
                  << (tri.vi[2] + 1) << "//" << ni << "\n";
            }
        }
    } else {
        f << "g default\n";
        for (std::size_t fi = 0; fi < m.faces.size(); ++fi) {
            const auto& tri = m.faces[fi];
            uint32_t ni = static_cast<uint32_t>(fi) + 1;
            f << "f "
              << (tri.vi[0] + 1) << "//" << ni << " "
              << (tri.vi[1] + 1) << "//" << ni << " "
              << (tri.vi[2] + 1) << "//" << ni << "\n";
        }
    }

    if (!f) throw std::runtime_error("OBJExporter: write error on " + path);
}

std::string OBJExporter::extension() const {
    return ".obj";
}

// ---------------------------------------------------------------------------
// PLYExporter
// ---------------------------------------------------------------------------

void PLYExporter::write(const SurfaceManifold& m, const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("PLYExporter: cannot open " + path);

    const bool has_color = (curvature_ != nullptr);
    const std::size_t nv = m.vertices.size();
    const std::size_t nf = m.faces.size();

    f << "ply\n";
    f << "format ascii 1.0\n";
    f << "element vertex " << nv << "\n";
    f << "property float x\n";
    f << "property float y\n";
    f << "property float z\n";
    if (has_color) {
        f << "property uchar red\n";
        f << "property uchar green\n";
        f << "property uchar blue\n";
    }
    f << "element face " << nf << "\n";
    f << "property list uchar int vertex_indices\n";
    f << "end_header\n";

    for (std::size_t vi = 0; vi < nv; ++vi) {
        const auto& pos = m.vertices[vi].pos;
        f << pos[0] << " " << pos[1] << " " << pos[2];
        if (has_color) {
            const auto& rgb = curvature_->jet_colors[vi];
            f << " " << static_cast<int>(rgb[0])
              << " " << static_cast<int>(rgb[1])
              << " " << static_cast<int>(rgb[2]);
        }
        f << "\n";
    }

    for (const auto& tri : m.faces) {
        f << "3 " << tri.vi[0] << " " << tri.vi[1] << " " << tri.vi[2] << "\n";
    }

    if (!f) throw std::runtime_error("PLYExporter: write error on " + path);
}

std::string PLYExporter::extension() const {
    return ".ply";
}

// ---------------------------------------------------------------------------
// STEPExporter
// NOTE: This is a simplified AP203 B-Rep stub.
// Production use requires Open CASCADE or equivalent STEP kernel.
// ---------------------------------------------------------------------------

void STEPExporter::write(const SurfaceManifold& m, const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("STEPExporter: cannot open " + path);

    f << "ISO-10303-21;\n";
    f << "HEADER;\n";
    f << "FILE_DESCRIPTION(('CFD Mesh Export'),'2;1');\n";
    f << "FILE_NAME('" << path << "','2026-05-10T00:00:00',(''),(''),'MeshMaster','','');\n";
    f << "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n";
    f << "ENDSEC;\n";
    f << "DATA;\n";

    uint32_t id = 1;

    for (const auto& tri : m.faces) {
        uint32_t pt_ids[3];
        for (int v = 0; v < 3; ++v) {
            const auto& pos = m.vertices[tri.vi[v]].pos;
            f << "#" << id << "=CARTESIAN_POINT('',("
              << pos[0] << "," << pos[1] << "," << pos[2] << "));\n";
            pt_ids[v] = id++;
        }

        uint32_t dir_id = id;
        f << "#" << id++ << "=DIRECTION('',("
          << tri.normal[0] << "," << tri.normal[1] << "," << tri.normal[2] << "));\n";

        f << "#" << id++ << "=ADVANCED_FACE('',(#" << pt_ids[0]
          << ",#" << pt_ids[1] << ",#" << pt_ids[2] << "),#" << dir_id << ",.T.);\n";
    }

    f << "ENDSEC;\n";
    f << "END-ISO-10303-21;\n";

    if (!f) throw std::runtime_error("STEPExporter: write error on " + path);
}

std::string STEPExporter::extension() const {
    return ".step";
}

// ---------------------------------------------------------------------------
// Factory
// ---------------------------------------------------------------------------

std::unique_ptr<MeshExporter> makeExporter(
    const std::string& format,
    const std::vector<std::vector<uint32_t>>* groups,
    const CurvatureResult*                    curvature)
{
    if (format == "stl_binary")
        return std::make_unique<STLExporter>(STLExporter::Mode::Binary);
    if (format == "stl_ascii")
        return std::make_unique<STLExporter>(STLExporter::Mode::ASCII);
    if (format == "obj")
        return std::make_unique<OBJExporter>(groups);
    if (format == "ply")
        return std::make_unique<PLYExporter>(curvature);
    if (format == "step")
        return std::make_unique<STEPExporter>();
    if (format == "vtu")
        return std::make_unique<VTUExporter>();

    throw std::invalid_argument("makeExporter: unknown format '" + format + "'");
}

// ---------------------------------------------------------------------------
// VTUExporter — VTK XML UnstructuredGrid
// ---------------------------------------------------------------------------

void VTUExporter::write(const SurfaceManifold& m, const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("VTUExporter: cannot open " + path);

    const size_t nv = m.vertices.size();
    const size_t nc = m.faces.size();

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << nv << "\" NumberOfCells=\"" << nc << "\">\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& v : m.vertices)
        f << "          " << v.pos[0] << " " << v.pos[1] << " " << v.pos[2] << "\n";
    f << "        </DataArray>\n      </Points>\n      <Cells>\n"
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& t : m.faces)
        f << "          " << t.vi[0] << " " << t.vi[1] << " " << t.vi[2] << "\n";
    f << "        </DataArray>\n"
      << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
    for (size_t i = 1; i <= nc; ++i) f << (i*3) << " ";
    f << "\n        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < nc; ++i) f << "5 ";  // VTK_TRIANGLE
    f << "\n        </DataArray>\n      </Cells>\n    </Piece>\n"
      << "  </UnstructuredGrid>\n</VTKFile>\n";
}

void VTUExporter::writeVolume(const VolumeMesh& vol, const std::string& path) const {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("VTUExporter: cannot open " + path);

    const size_t nv = vol.vertices.size();
    const size_t nt = vol.tets.size();
    const size_t np = vol.prisms.size();
    const size_t nc = nt + np;

    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << nv << "\" NumberOfCells=\"" << nc << "\">\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& v : vol.vertices)
        f << "          " << v.pos[0] << " " << v.pos[1] << " " << v.pos[2] << "\n";
    f << "        </DataArray>\n      </Points>\n";

    // Cell data: cell_type and layer
    f << "      <CellData>\n"
      << "        <DataArray type=\"Int32\" Name=\"cell_type\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < nt; ++i) f << "0 ";
    for (size_t i = 0; i < np; ++i) f << "1 ";
    f << "\n        </DataArray>\n"
      << "        <DataArray type=\"Int32\" Name=\"layer\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < nt; ++i) f << "-1 ";
    for (const auto& pr : vol.prisms) f << pr.layer << " ";
    f << "\n        </DataArray>\n      </CellData>\n";

    // Connectivity
    f << "      <Cells>\n"
      << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& t : vol.tets)
        f << "          " << t.nodes[0] << " " << t.nodes[1] << " " << t.nodes[2] << " " << t.nodes[3] << "\n";
    for (const auto& pr : vol.prisms)
        f << "          " << pr.nodes[0] << " " << pr.nodes[1] << " " << pr.nodes[2]
          << " " << pr.nodes[3] << " " << pr.nodes[4] << " " << pr.nodes[5] << "\n";
    f << "        </DataArray>\n";

    // Offsets
    f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ";
    int off = 0;
    for (size_t i = 0; i < nt; ++i){ off += 4; f << off << " "; }
    for (size_t i = 0; i < np; ++i){ off += 6; f << off << " "; }
    f << "\n        </DataArray>\n";

    // Types: 10=VTK_TETRA, 13=VTK_WEDGE
    f << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
    for (size_t i = 0; i < nt; ++i) f << "10 ";
    for (size_t i = 0; i < np; ++i) f << "13 ";
    f << "\n        </DataArray>\n      </Cells>\n    </Piece>\n"
      << "  </UnstructuredGrid>\n</VTKFile>\n";
}

} // namespace cfd
