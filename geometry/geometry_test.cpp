//
// geometry_test.cpp — self-contained unit tests for GeometryEngine
// Compile: see CMakeLists.txt
//
#include "geometry_engine.hpp"

#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ── helpers ─────────────────────────────────────────────────────────────────

static void CHECK(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        std::exit(1);
    }
}

static bool approx_eq(double a, double b, double eps = 1e-6) {
    return std::abs(a - b) < eps;
}

// Write a minimal binary STL containing `tris` triangles.
// Layout: 80-byte header | uint32 count | N × (3f normal + 9f verts + u16 attr)
static void writeBinarySTL(const std::string& path,
                            const std::vector<cfd::Triangle>& tris,
                            const std::vector<cfd::Vertex>&   verts) {
    std::ofstream f(path, std::ios::binary);

    // 80-byte header
    char header[80] = {};
    std::strncpy(header, "Binary STL test", sizeof(header));
    f.write(header, 80);

    auto count = static_cast<uint32_t>(tris.size());
    f.write(reinterpret_cast<const char*>(&count), 4);

    for (const auto& tri : tris) {
        auto wf = [&](double d) {
            float fv = static_cast<float>(d);
            f.write(reinterpret_cast<const char*>(&fv), 4);
        };
        wf(tri.normal[0]); wf(tri.normal[1]); wf(tri.normal[2]);
        for (int v = 0; v < 3; ++v) {
            const auto& p = verts[tri.vi[v]].pos;
            wf(p[0]); wf(p[1]); wf(p[2]);
        }
        uint16_t attr = 0;
        f.write(reinterpret_cast<const char*>(&attr), 2);
    }
}

// Write an ASCII STL file
static void writeAsciiSTL(const std::string& path,
                           const std::vector<cfd::Triangle>& tris,
                           const std::vector<cfd::Vertex>&   verts) {
    std::ofstream f(path);
    f << "solid test\n";
    for (const auto& tri : tris) {
        f << "  facet normal "
          << tri.normal[0] << " "
          << tri.normal[1] << " "
          << tri.normal[2] << "\n";
        f << "    outer loop\n";
        for (int v = 0; v < 3; ++v) {
            const auto& p = verts[tri.vi[v]].pos;
            f << "      vertex " << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        f << "    endloop\n";
        f << "  endfacet\n";
    }
    f << "endsolid test\n";
}

// Build a closed tetrahedron mesh (4 triangles, 4 vertices, watertight)
static std::pair<std::vector<cfd::Triangle>, std::vector<cfd::Vertex>>
makeTetrahedron() {
    using namespace cfd;
    std::vector<Vertex> V = {
        Vertex{{ 1.0,  1.0,  1.0}},
        Vertex{{-1.0, -1.0,  1.0}},
        Vertex{{-1.0,  1.0, -1.0}},
        Vertex{{ 1.0, -1.0, -1.0}}
    };

    // Four faces of the tetrahedron (consistent outward winding)
    std::vector<Triangle> T = {
        Triangle{{ {0, 1, 2} }, {0,0,0}},
        Triangle{{ {0, 3, 1} }, {0,0,0}},
        Triangle{{ {0, 2, 3} }, {0,0,0}},
        Triangle{{ {1, 3, 2} }, {0,0,0}}
    };
    return {T, V};
}

// ── tests ────────────────────────────────────────────────────────────────────

void test_binary_stl_load(const std::string& tmp_dir) {
    auto [tris, verts] = makeTetrahedron();
    std::string path = tmp_dir + "/tet.stl";
    writeBinarySTL(path, tris, verts);

    cfd::GeometryEngine eng;
    auto manifold = eng.importCAD(path);

    CHECK(manifold.faces.size() == 4,    "binary STL: 4 faces");
    CHECK(manifold.vertices.size() == 12,"binary STL: 12 raw vertices (pre-merge)");
    std::cout << "PASS: binary_stl_load\n";
}

void test_ascii_stl_load(const std::string& tmp_dir) {
    auto [tris, verts] = makeTetrahedron();
    std::string path = tmp_dir + "/tet_ascii.stl";
    writeAsciiSTL(path, tris, verts);

    cfd::GeometryEngine eng;
    auto manifold = eng.importCAD(path);

    CHECK(manifold.faces.size() == 4,  "ascii STL: 4 faces");
    CHECK(!manifold.faces.empty(),     "ascii STL: non-empty");
    std::cout << "PASS: ascii_stl_load\n";
}

void test_merge_vertices(const std::string& tmp_dir) {
    auto [tris, verts] = makeTetrahedron();
    std::string path = tmp_dir + "/tet_merge.stl";
    writeBinarySTL(path, tris, verts);

    cfd::GeometryEngine eng;
    auto manifold = eng.importCAD(path);

    // Before repair: 12 vertices (3 per face, shared vertices duplicated)
    CHECK(manifold.vertices.size() == 12, "pre-merge: 12 verts");

    bool ok = eng.repairManifold(manifold);

    // After merge: 4 unique vertices
    CHECK(manifold.vertices.size() == 4, "post-merge: 4 unique verts");
    CHECK(ok, "tetrahedron is watertight after repair");
    std::cout << "PASS: merge_vertices\n";
}

void test_watertight_check(const std::string& tmp_dir) {
    auto [tris, verts] = makeTetrahedron();
    std::string path = tmp_dir + "/tet_wt.stl";
    writeBinarySTL(path, tris, verts);

    cfd::GeometryEngine eng;
    auto manifold = eng.importCAD(path);
    eng.repairManifold(manifold);

    CHECK(manifold.is_watertight, "tetrahedron is watertight");
    std::cout << "PASS: watertight_check\n";
}

void test_fill_holes(const std::string& tmp_dir) {
    // Build a tetrahedron but omit one face → open hole
    auto [tris, verts] = makeTetrahedron();
    std::vector<cfd::Triangle> open_tris = {tris[0], tris[1], tris[2]};  // drop face 3
    std::string path = tmp_dir + "/open_tet.stl";
    writeBinarySTL(path, open_tris, verts);

    cfd::GeometryEngine eng;
    auto manifold = eng.importCAD(path);

    CHECK(manifold.faces.size() == 3, "open tet: 3 faces loaded");

    bool ok = eng.repairManifold(manifold);

    // fillHoles should patch the 3-edge boundary loop
    CHECK(manifold.faces.size() > 3, "fill_holes: new triangles added");
    CHECK(ok, "patched manifold is watertight");
    std::cout << "PASS: fill_holes\n";
}

void test_identify_boundaries_default(const std::string& tmp_dir) {
    auto [tris, verts] = makeTetrahedron();
    std::string path = tmp_dir + "/tet_bc.stl";
    writeBinarySTL(path, tris, verts);

    cfd::GeometryEngine eng;
    auto manifold = eng.importCAD(path);
    eng.repairManifold(manifold);

    // No named patches → identifyBoundaries should create one "wall" patch
    eng.identifyBoundaries(manifold, {});

    CHECK(manifold.patches.size() == 1, "one default patch");
    CHECK(manifold.patches[0].name == "wall", "default patch named 'wall'");
    CHECK(manifold.patches[0].type == cfd::BoundaryType::Wall,
          "default patch type Wall");
    std::cout << "PASS: identify_boundaries_default\n";
}

void test_identify_boundaries_rules() {
    // Manually build a manifold with named patches and check rule matching
    cfd::SurfaceManifold m;
    cfd::BoundaryPatch p1, p2, p3;
    p1.name = "inlet_top";
    p2.name = "outlet_face";
    p3.name = "symmetry_plane";

    m.patches = {p1, p2, p3};

    cfd::GeometryEngine eng;
    eng.identifyBoundaries(m, {});  // use built-in heuristics

    CHECK(m.patches[0].type == cfd::BoundaryType::Inlet,    "inlet_top → Inlet");
    CHECK(m.patches[1].type == cfd::BoundaryType::Outlet,   "outlet_face → Outlet");
    CHECK(m.patches[2].type == cfd::BoundaryType::Symmetry, "symmetry_plane → Symmetry");
    std::cout << "PASS: identify_boundaries_rules\n";
}

void test_patch_centroid(const std::string& tmp_dir) {
    // Single equilateral triangle in XY plane: centroid should be at (1/3, 1/3, 0)
    cfd::SurfaceManifold m;
    m.vertices = {
        cfd::Vertex{{0.0, 0.0, 0.0}},
        cfd::Vertex{{1.0, 0.0, 0.0}},
        cfd::Vertex{{0.0, 1.0, 0.0}}
    };
    cfd::Triangle tri;
    tri.vi     = {0, 1, 2};
    tri.normal = {0.0, 0.0, 1.0};
    m.faces.push_back(tri);

    cfd::BoundaryPatch patch;
    patch.name = "test";
    patch.type = cfd::BoundaryType::Wall;
    patch.face_ids = {0};
    m.patches.push_back(patch);

    cfd::GeometryEngine eng;
    auto c = eng.patchCentroid(m, patch);

    CHECK(approx_eq(c[0], 1.0/3.0), "centroid x ≈ 1/3");
    CHECK(approx_eq(c[1], 1.0/3.0), "centroid y ≈ 1/3");
    CHECK(approx_eq(c[2], 0.0),     "centroid z = 0");
    std::cout << "PASS: patch_centroid\n";
}

void test_progress_callback(const std::string& tmp_dir) {
    auto [tris, verts] = makeTetrahedron();
    std::string path = tmp_dir + "/tet_prog.stl";
    writeBinarySTL(path, tris, verts);

    std::vector<std::pair<std::string, double>> log;
    cfd::GeometryEngine eng([&](const std::string& stage, double frac) {
        log.push_back({stage, frac});
    });

    auto manifold = eng.importCAD(path);
    eng.repairManifold(manifold);

    CHECK(!log.empty(), "progress callback was invoked");
    // Last import event should be at fraction 1.0
    auto it = std::find_if(log.begin(), log.end(),
                           [](const auto& e){ return e.first == "import:done"; });
    CHECK(it != log.end(), "import:done event fired");
    CHECK(approx_eq(it->second, 1.0), "import:done fraction == 1.0");
    std::cout << "PASS: progress_callback\n";
}

void test_bad_file_throws() {
    cfd::GeometryEngine eng;
    bool threw = false;
    try {
        eng.importCAD("/tmp/__nonexistent_file_xyz__.stl");
    } catch (const cfd::GeometryError&) {
        threw = true;
    }
    CHECK(threw, "importCAD throws GeometryError for missing file");
    std::cout << "PASS: bad_file_throws\n";
}

// ── main ─────────────────────────────────────────────────────────────────────

int main() {
    const std::string tmp = "/tmp/cfd_geometry_tests";
    // Create tmp directory
    std::system(("mkdir -p " + tmp).c_str());

    std::cout << "=== GeometryEngine unit tests ===\n";

    test_binary_stl_load(tmp);
    test_ascii_stl_load(tmp);
    test_merge_vertices(tmp);
    test_watertight_check(tmp);
    test_fill_holes(tmp);
    test_identify_boundaries_default(tmp);
    test_identify_boundaries_rules();
    test_patch_centroid(tmp);
    test_progress_callback(tmp);
    test_bad_file_throws();

    std::cout << "\nAll tests passed.\n";
    return 0;
}
