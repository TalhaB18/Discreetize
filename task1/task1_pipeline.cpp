// Task 1 Pipeline — Boolean Subtraction + Tetrahedral Meshing + All-Quad Surface Mesh
// Intern Interview, Discreetize Computational Geometry
//
// Usage:
//   task1_pipeline [--output-dir <dir>] [--resolution <int>] [--no-download]
//
// Downloads 5 Thingi10K models from HuggingFace, performs boolean subtraction
// for 3 representative pairs, generates tetrahedral volume meshes, and exports
// all-quad surface meshes (bonus).

#include "../geometry/geometry_engine.hpp"
#include "../geometry/mesh_sdf.hpp"
#include "../geometry/primitives.hpp"
#include "../meshing/mesh_engine.hpp"
#include "../meshing/quad_mesh.hpp"
#include "../export/mesh_exporter.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <string>
#include <vector>
#include <sys/stat.h>

using namespace cfd;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static void ensure_dir(const std::string& d){
    std::string cmd = "mkdir -p \"" + d + "\"";
    std::system(cmd.c_str());
}

static bool file_exists(const std::string& p){
    struct stat s{};
    return stat(p.c_str(), &s) == 0;
}

static void download_model(const std::string& url, const std::string& dest){
    if(file_exists(dest)){
        std::printf("  [cached] %s\n", dest.c_str());
        return;
    }
    std::printf("  Downloading %s ...\n", url.c_str());
    std::string cmd = "curl -L -f --silent --show-error -o \"" + dest + "\" \"" + url + "\"";
    int rc = std::system(cmd.c_str());
    if(rc != 0) std::fprintf(stderr, "  WARNING: download failed for %s\n", url.c_str());
    else std::printf("  Saved: %s\n", dest.c_str());
}

// Normalize a surface mesh: translate centroid to origin, scale to fit [-1,1]^3
static SurfaceManifold normalize_mesh(SurfaceManifold m){
    if(m.vertices.empty()) return m;
    double cx=0,cy=0,cz=0;
    for(const auto& v : m.vertices){ cx+=v.pos[0]; cy+=v.pos[1]; cz+=v.pos[2]; }
    cx/=m.vertices.size(); cy/=m.vertices.size(); cz/=m.vertices.size();
    double r=0;
    for(auto& v : m.vertices){
        v.pos[0]-=cx; v.pos[1]-=cy; v.pos[2]-=cz;
        double d=std::sqrt(v.pos[0]*v.pos[0]+v.pos[1]*v.pos[1]+v.pos[2]*v.pos[2]);
        if(d>r) r=d;
    }
    if(r>1e-12){ for(auto& v : m.vertices){ v.pos[0]/=r; v.pos[1]/=r; v.pos[2]/=r; } }
    // Update normals
    for(auto& tri : m.faces){
        auto& A=m.vertices[tri.vi[0]].pos;
        auto& B=m.vertices[tri.vi[1]].pos;
        auto& C=m.vertices[tri.vi[2]].pos;
        double ux=B[0]-A[0],uy=B[1]-A[1],uz=B[2]-A[2];
        double vx=C[0]-A[0],vy=C[1]-A[1],vz=C[2]-A[2];
        double nx=uy*vz-uz*vy,ny=uz*vx-ux*vz,nz=ux*vy-uy*vx;
        double len=std::sqrt(nx*nx+ny*ny+nz*nz);
        if(len>1e-12){ tri.normal={nx/len,ny/len,nz/len}; }
    }
    return m;
}

// Scale and shift mesh B for a good boolean subtraction overlap
static SurfaceManifold transform_mesh(SurfaceManifold m, double scale, double tx, double ty, double tz){
    for(auto& v : m.vertices){
        v.pos[0]=v.pos[0]*scale+tx;
        v.pos[1]=v.pos[1]*scale+ty;
        v.pos[2]=v.pos[2]*scale+tz;
    }
    return m;
}

// Boolean subtract A - B using SDF + Marching Cubes
static SurfaceManifold boolean_subtract_meshes(
    const SurfaceManifold& A, const SurfaceManifold& B, int resolution)
{
    MeshSDF sdf_a(A), sdf_b(B);

    Vec3d lo_a, hi_a, lo_b, hi_b;
    sdf_a.bounds(lo_a, hi_a);
    sdf_b.bounds(lo_b, hi_b);

    double lo[3] = {std::min(lo_a[0],lo_b[0]), std::min(lo_a[1],lo_b[1]), std::min(lo_a[2],lo_b[2])};
    double hi[3] = {std::max(hi_a[0],hi_b[0]), std::max(hi_a[1],hi_b[1]), std::max(hi_a[2],hi_b[2])};

    auto sdf = [&](double x, double y, double z) -> double {
        return std::max(sdf_a.eval({x,y,z}), -sdf_b.eval({x,y,z}));
    };

    std::printf("    Running Marching Cubes (resolution %d)...\n", resolution);
    auto tris = marching_cubes_run(sdf, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], resolution);
    std::printf("    MC triangles: %zu\n", tris.size()/3);

    auto result = tris_to_manifold(tris);
    GeometryEngine geo;
    geo.mergeVertices(result, 1e-8);
    geo.repairManifold(result);
    return result;
}

// Process one pair: subtract, tet mesh, quad mesh, export
static void process_pair(
    const std::string& label,
    SurfaceManifold mesh_a,
    SurfaceManifold mesh_b,
    const std::string& out_dir,
    int mc_resolution)
{
    std::printf("\n--- %s ---\n", label.c_str());

    auto t0 = std::chrono::steady_clock::now();

    // Boolean subtraction
    std::printf("  Boolean subtraction (A - B)...\n");
    auto result = boolean_subtract_meshes(mesh_a, mesh_b, mc_resolution);
    std::printf("  Result: %zu vertices, %zu faces\n",
                result.vertices.size(), result.faces.size());

    auto t1 = std::chrono::steady_clock::now();
    double dt_bool = std::chrono::duration<double>(t1-t0).count();
    std::printf("  Boolean time: %.2f s\n", dt_bool);

    if(result.faces.empty()){
        std::printf("  WARNING: empty result, skipping meshing.\n");
        return;
    }

    // Export boolean result as STL + VTU
    std::string stl_path = out_dir + "/" + label + "_boolean.stl";
    std::string vtu_surf = out_dir + "/" + label + "_boolean.vtu";
    makeExporter("stl_binary")->write(result, stl_path);
    VTUExporter vtu_exp;
    vtu_exp.write(result, vtu_surf);
    std::printf("  Saved: %s\n", stl_path.c_str());
    std::printf("  Saved: %s\n", vtu_surf.c_str());

    // Tag all faces as Wall for BL inflation
    {
        BoundaryPatch wall; wall.name="wall"; wall.type=BoundaryType::Wall;
        wall.face_ids.resize(result.faces.size());
        for(uint32_t i=0;i<(uint32_t)result.faces.size();++i) wall.face_ids[i]=i;
        result.patches.push_back(std::move(wall));
    }

    // Tetrahedral mesh
    std::printf("  Generating tetrahedral mesh...\n");
    MeshEngine eng([](const std::string& s, double p){
        std::printf("    [%3.0f%%] %s\n", p*100.0, s.c_str());
    });
    BLParameters bl;
    bl.first_cell_height = 0.02;
    bl.growth_ratio      = 1.2;
    bl.num_layers        = 5;
    double rn = std::pow(bl.growth_ratio, bl.num_layers);
    bl.total_thickness = bl.first_cell_height * (rn - 1.0) / (bl.growth_ratio - 1.0);

    try {
        auto prisms = eng.inflatePrisms(result, bl);
        auto vol    = eng.generateTetCore(prisms, 0.08);
        auto qual   = eng.checkQuality(vol);

        std::string vtu_vol = out_dir + "/" + label + "_tet.vtu";
        vtu_exp.writeVolume(vol, vtu_vol);
        std::printf("  Saved: %s\n", vtu_vol.c_str());

        std::printf("  Tet mesh — %zu tets, %zu prisms, %zu vertices\n",
                    vol.tets.size(), vol.prisms.size(), vol.vertices.size());
        std::printf("  Quality — max_skewness: %.3f  avg_skewness: %.3f  bad_cells: %d\n",
                    qual.max_skewness, qual.avg_skewness, qual.bad_cells);
    } catch(const std::exception& e){
        std::fprintf(stderr, "  Meshing error: %s\n", e.what());
    }

    auto t2 = std::chrono::steady_clock::now();
    double dt_mesh = std::chrono::duration<double>(t2-t1).count();
    std::printf("  Mesh time: %.2f s\n", dt_mesh);

    // All-quad surface mesh (bonus)
    std::printf("  Generating all-quad surface mesh (bonus)...\n");
    QuadMesher qm;
    auto quad = qm.convert(result);
    std::string quad_path = out_dir + "/" + label + "_quad.vtu";
    qm.writeVTU(quad, quad_path);
    std::printf("  Saved: %s  (%zu quads)\n", quad_path.c_str(), quad.quads.size());

    double dt_total = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::printf("  Total: %.2f s\n", dt_total);
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    std::string out_dir = "output/task1";
    int mc_resolution   = 40;
    bool no_download    = false;

    for(int i=1; i<argc; ++i){
        if(!std::strcmp(argv[i],"--output-dir") && i+1<argc) out_dir = argv[++i];
        else if(!std::strcmp(argv[i],"--resolution") && i+1<argc) mc_resolution = std::atoi(argv[++i]);
        else if(!std::strcmp(argv[i],"--no-download")) no_download = true;
        else if(!std::strcmp(argv[i],"--help")){
            std::printf("Usage: task1_pipeline [--output-dir <dir>] [--resolution <int>] [--no-download]\n");
            return 0;
        }
    }

    std::printf("=== Task 1: Boolean Subtraction + Tetrahedral Meshing ===\n");
    std::printf("Output dir:  %s\n", out_dir.c_str());
    std::printf("MC resolution: %d\n\n", mc_resolution);

    std::string model_dir = out_dir + "/models";
    ensure_dir(model_dir);

    // Model URLs and local paths
    struct Model { std::string id, url, path; };
    const std::string BASE = "https://huggingface.co/datasets/Thingi10K/Thingi10K/resolve/main/raw_meshes/";
    std::vector<Model> models = {
        {"100035",  BASE+"100035.stl",  model_dir+"/100035.stl"},
        {"100075",  BASE+"100075.stl",  model_dir+"/100075.stl"},
        {"100335",  BASE+"100335.stl",  model_dir+"/100335.stl"},
        {"100388",  BASE+"100388.stl",  model_dir+"/100388.stl"},
        {"1005587", BASE+"1005587.stl", model_dir+"/1005587.stl"},
    };

    if(!no_download){
        std::printf("--- Downloading Thingi10K models ---\n");
        for(auto& m : models) download_model(m.url, m.path);
    }

    // Load and repair all models
    std::printf("\n--- Loading and repairing models ---\n");
    GeometryEngine geo;
    std::vector<SurfaceManifold> meshes(models.size());
    for(size_t i=0; i<models.size(); ++i){
        if(!file_exists(models[i].path)){
            std::fprintf(stderr, "Model not found: %s\n", models[i].path.c_str());
            continue;
        }
        std::printf("  Loading %s...\n", models[i].id.c_str());
        try {
            meshes[i] = geo.importCAD(models[i].path);
            geo.repairManifold(meshes[i]);
            meshes[i] = normalize_mesh(std::move(meshes[i]));
            std::printf("    %zu vertices, %zu faces, watertight: %s\n",
                        meshes[i].vertices.size(), meshes[i].faces.size(),
                        meshes[i].is_watertight ? "yes" : "no");
        } catch(const std::exception& e){
            std::fprintf(stderr, "    Error: %s\n", e.what());
        }
    }

    // Pair 1: 100035 (Tail Section) minus 100075 (Robot Arm) scaled to 60%
    if(!meshes[0].faces.empty() && !meshes[1].faces.empty()){
        auto B = transform_mesh(meshes[1], 0.60, 0.30, 0.30, 0.30);
        process_pair("pair1_tail_minus_arm", meshes[0], B, out_dir, mc_resolution);
    }

    // Pair 2: 100335 (Stonehenge) minus 100388 (Menou) scaled to 50%
    if(!meshes[2].faces.empty() && !meshes[3].faces.empty()){
        auto B = transform_mesh(meshes[3], 0.50, 0.20, 0.00, 0.20);
        process_pair("pair2_stonehenge_minus_menou", meshes[2], B, out_dir, mc_resolution);
    }

    // Pair 3: 1005587 (Bird Feeder) minus 100035 (Tail Section) scaled to 40%
    if(!meshes[4].faces.empty() && !meshes[0].faces.empty()){
        auto B = transform_mesh(meshes[0], 0.40, 0.00, 0.00, 0.00);
        process_pair("pair3_birdfeeder_minus_tail", meshes[4], B, out_dir, mc_resolution);
    }

    std::printf("\n=== Task 1 complete. Results in: %s ===\n", out_dir.c_str());
    return 0;
}
