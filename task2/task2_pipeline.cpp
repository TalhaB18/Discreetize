// Task 2 Pipeline — Twisted Bifilar Solid: Generate + BL Mesh + Export
// Intern Interview, Discreetize Computational Geometry
//
// Usage:
//   task2_pipeline [--R <radius>] [--H <height>] [--revs <n>]
//                  [--resolution <int>] [--output-dir <dir>]
//
// Generates the twisted bifilar fluid domain (two helical sub-cylinders inside
// an outer cylinder), meshes it with prismatic boundary layer inflation and a
// tetrahedral core, and exports VTU files for ParaView visualisation.

#include "../geometry/bifilar.hpp"
#include "../geometry/geometry_engine.hpp"
#include "../meshing/mesh_engine.hpp"
#include "../export/mesh_exporter.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <string>
#include <sys/stat.h>

using namespace cfd;

static void ensure_dir(const std::string& d){
    std::string cmd = "mkdir -p \"" + d + "\"";
    std::system(cmd.c_str());
}

int main(int argc, char* argv[])
{
    BifilarConfig cfg;
    std::string out_dir = "output/task2";

    for(int i=1; i<argc; ++i){
        if(!std::strcmp(argv[i],"--R") && i+1<argc)          cfg.R          = std::atof(argv[++i]);
        else if(!std::strcmp(argv[i],"--H") && i+1<argc)     cfg.H          = std::atof(argv[++i]);
        else if(!std::strcmp(argv[i],"--revs") && i+1<argc)  cfg.revolutions= std::atof(argv[++i]);
        else if(!std::strcmp(argv[i],"--resolution") && i+1<argc) cfg.resolution = std::atoi(argv[++i]);
        else if(!std::strcmp(argv[i],"--output-dir") && i+1<argc) out_dir    = argv[++i];
        else if(!std::strcmp(argv[i],"--help")){
            std::printf("Usage: task2_pipeline [--R <r>] [--H <h>] [--revs <n>] "
                        "[--resolution <int>] [--output-dir <dir>]\n");
            return 0;
        }
    }

    ensure_dir(out_dir);

    std::printf("=== Task 2: Twisted Bifilar Solid ===\n");
    std::printf("  R = %.3f  r = %.3f  offset = %.3f\n", cfg.R, cfg.r, cfg.offset);
    std::printf("  H = %.3f  revolutions = %.1f\n", cfg.H, cfg.revolutions);
    std::printf("  MC resolution = %d  output = %s\n\n", cfg.resolution, out_dir.c_str());

    // ── 1. Generate surface ───────────────────────────────────────────────────
    std::printf("Step 1: Generating bifilar surface via Marching Cubes...\n");
    auto t0 = std::chrono::steady_clock::now();

    SurfaceManifold surface = generateBifilarSurface(cfg);

    double dt_gen = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::printf("  Vertices: %zu   Faces: %zu   (%.2f s)\n",
                surface.vertices.size(), surface.faces.size(), dt_gen);
    std::printf("  Watertight: %s\n", surface.is_watertight ? "YES" : "NO");

    // ── 2. Export surface STL ─────────────────────────────────────────────────
    std::string stl_path = out_dir + "/bifilar_surface.stl";
    makeExporter("stl_binary")->write(surface, stl_path);
    std::printf("  Saved surface: %s\n", stl_path.c_str());

    // ── 3. Export surface VTU (for ParaView) ──────────────────────────────────
    VTUExporter vtu;
    std::string surf_vtu = out_dir + "/bifilar_surface.vtu";
    vtu.write(surface, surf_vtu);
    std::printf("  Saved surface VTU: %s\n\n", surf_vtu.c_str());

    // Tag all faces as Wall so inflatePrisms can find them
    {
        BoundaryPatch wall;
        wall.name = "wall";
        wall.type = BoundaryType::Wall;
        wall.face_ids.resize(surface.faces.size());
        for(uint32_t i=0; i<(uint32_t)surface.faces.size(); ++i) wall.face_ids[i]=i;
        surface.patches.push_back(std::move(wall));
    }

    // ── 4. Boundary layer prism inflation ────────────────────────────────────
    std::printf("Step 2: Prismatic boundary layer inflation...\n");
    BLParameters bl;
    bl.first_cell_height = 0.008 * cfg.R;
    bl.growth_ratio      = 1.1;
    bl.num_layers        = 10;
    double rn = std::pow(bl.growth_ratio, bl.num_layers);
    bl.total_thickness = bl.first_cell_height * (rn - 1.0) / (bl.growth_ratio - 1.0);

    std::printf("  First layer height δ = %.5f R\n", bl.first_cell_height / cfg.R);
    std::printf("  Layers: %d   Growth ratio: %.2f\n", bl.num_layers, bl.growth_ratio);
    std::printf("  Total BL thickness: %.5f\n", bl.total_thickness);

    MeshEngine eng([](const std::string& s, double p){
        std::printf("  [%3.0f%%] %s\n", p*100.0, s.c_str());
    });

    auto t1 = std::chrono::steady_clock::now();
    PrismMesh prisms;
    try {
        prisms = eng.inflatePrisms(surface, bl);
    } catch(const std::exception& e){
        std::fprintf(stderr, "BL inflation error: %s\n", e.what());
        return 1;
    }
    double dt_bl = std::chrono::duration<double>(std::chrono::steady_clock::now()-t1).count();
    std::printf("  Prism elements: %zu   (%.2f s)\n\n", prisms.prisms.size(), dt_bl);

    // ── 5. Tetrahedral core ───────────────────────────────────────────────────
    std::printf("Step 3: Generating tetrahedral core...\n");
    double core_size = 0.06 * cfg.R;
    std::printf("  Target cell size: %.4f\n", core_size);

    auto t2 = std::chrono::steady_clock::now();
    VolumeMesh vol;
    try {
        vol = eng.generateTetCore(prisms, core_size);
    } catch(const std::exception& e){
        std::fprintf(stderr, "Tet core error: %s\n", e.what());
        return 1;
    }
    double dt_tet = std::chrono::duration<double>(std::chrono::steady_clock::now()-t2).count();
    std::printf("  Tets: %zu   Prisms: %zu   Vertices: %zu   (%.2f s)\n\n",
                vol.tets.size(), vol.prisms.size(), vol.vertices.size(), dt_tet);

    // ── 6. Quality report ────────────────────────────────────────────────────
    std::printf("Step 4: Mesh quality assessment...\n");
    auto qual = eng.checkQuality(vol);
    std::printf("  Max skewness:     %.4f  (threshold 0.85)\n", qual.max_skewness);
    std::printf("  Avg skewness:     %.4f\n",                    qual.avg_skewness);
    std::printf("  Max aspect ratio: %.4f\n",                    qual.max_aspect_ratio);
    std::printf("  Min orthogonality:%.2f deg\n",                qual.min_orthogonality);
    std::printf("  Bad cells:        %d\n\n",                    qual.bad_cells);

    // ── 7. Export volume VTU ─────────────────────────────────────────────────
    std::printf("Step 5: Exporting mesh...\n");
    std::string vol_vtu = out_dir + "/bifilar_tet_mesh.vtu";
    vtu.writeVolume(vol, vol_vtu);
    std::printf("  Saved volume VTU: %s\n", vol_vtu.c_str());

    double dt_total = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::printf("\n=== Task 2 complete in %.2f s. Results in: %s ===\n", dt_total, out_dir.c_str());
    return 0;
}
