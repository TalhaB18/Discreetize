// MeshMaster v2.0 — 3D Mesh Processing Tool
#include "geometry/geometry_engine.hpp"
#include "analysis/curvature.hpp"
#include "slicer/slicer.hpp"
#include "export/mesh_exporter.hpp"
#include <iostream>
#include <string>
#include <algorithm>

// ---------------------------------------------------------------------------
// Theme — ANSI colour codes
// ---------------------------------------------------------------------------
struct Theme {
    const char* reset;
    const char* title;
    const char* ok;
    const char* warn;
    const char* err;
    const char* label;
    const char* dim;
};

static const Theme DARK = {
    "\033[0m",
    "\033[1;36m",   // bold cyan
    "\033[1;32m",   // bold green
    "\033[1;33m",   // bold yellow
    "\033[1;31m",   // bold red
    "\033[1;37m",   // bold white
    "\033[2;37m",   // dim grey
};

static const Theme LIGHT = {
    "\033[0m",
    "\033[1;34m",   // bold blue
    "\033[0;32m",   // green
    "\033[0;33m",   // dark yellow
    "\033[0;31m",   // red
    "\033[1;30m",   // bold dark grey
    "\033[2;30m",   // dim black
};

static const Theme* T = &DARK;

// ---------------------------------------------------------------------------
// Banner
// ---------------------------------------------------------------------------
static void banner() {
    std::cout
        << T->title
        << "\n"
           "  ███╗   ███╗███████╗███████╗██╗  ██╗\n"
           "  ████╗ ████║██╔════╝██╔════╝██║  ██║\n"
           "  ██╔████╔██║█████╗  ███████╗███████║\n"
           "  ██║╚██╔╝██║██╔══╝  ╚════██║██╔══██║\n"
           "  ██║ ╚═╝ ██║███████╗███████║██║  ██║\n"
           "  ╚═╝     ╚═╝╚══════╝╚══════╝╚═╝  ╚═╝\n"
        << T->dim  << "  MeshMaster v2.0  |  3D Mesh Processing Tool\n"
        << T->reset << "\n";
}

// ---------------------------------------------------------------------------
// Help
// ---------------------------------------------------------------------------
static void printHelp() {
    banner();
    std::cout
        << T->label << "Usage: " << T->reset
        << "meshmaster [--dark|--light] <command> <input.stl> [options]\n\n"
        << T->title << "Commands:\n" << T->reset
        << "  " << T->label << "info"      << T->reset << "      <input.stl>\n"
        << "               Vertex/face count, watertight check, manifold report.\n\n"
        << "  " << T->label << "heal"      << T->reset << "      <input.stl> [-o output.stl]\n"
        << "               Merge duplicate verts, fix normals, fill holes. Default: healed.stl\n\n"
        << "  " << T->label << "slice"     << T->reset << "      <input.stl> -z <value> [-o prefix]\n"
        << "               Slice at Z → saves <prefix>.svg and <prefix>.txt\n\n"
        << "  " << T->label << "curvature" << T->reset << " <input.stl> [-o output.ply]\n"
        << "               Mean curvature heatmap → PLY with per-vertex colour.\n\n"
        << "  " << T->label << "export"    << T->reset << "    <input.stl> -f <format> [-o output]\n"
        << "               Formats: " << T->dim << "stl  stl_ascii  obj  ply  step" << T->reset << "\n\n"
        << T->title << "Theme flags:\n" << T->reset
        << "  --dark    Dark terminal background (default)\n"
        << "  --light   Light terminal background\n\n"
        << T->title << "Examples:\n" << T->reset
        << T->dim
        << "  meshmaster heal      my_model.stl -o fixed.stl\n"
        << "  meshmaster slice     my_model.stl -z 12.5 -o crosssection\n"
        << "  meshmaster curvature my_model.stl -o heat.ply\n"
        << "  meshmaster export    my_model.stl -f obj -o model.obj\n"
        << T->reset;
}

// ---------------------------------------------------------------------------
// Argument helpers
// ---------------------------------------------------------------------------
static std::string getArg(int argc, char* argv[], const std::string& flag, const std::string& def = "") {
    for (int i = 1; i < argc - 1; ++i)
        if (argv[i] == flag) return argv[i + 1];
    return def;
}

static bool hasFlag(int argc, char* argv[], const std::string& flag) {
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == flag) return true;
    return false;
}

// ---------------------------------------------------------------------------
// Load + heal helper
// ---------------------------------------------------------------------------
static cfd::SurfaceManifold loadAndHeal(const std::string& path, bool quiet = false) {
    auto progress = [&](const std::string& stage, double pct) {
        if (!quiet)
            std::cout << T->dim << "  [" << int(pct * 100) << "%] " << stage << T->reset << "\n";
    };
    cfd::GeometryEngine geo(progress);
    auto surface = geo.importCAD(path);
    surface.is_watertight = geo.repairManifold(surface);
    return surface;
}

// ---------------------------------------------------------------------------
// info
// ---------------------------------------------------------------------------
static int cmd_info(const std::string& path) {
    std::cout << T->label << "Loading: " << T->reset << path << "\n";
    cfd::GeometryEngine geo;
    auto surface = geo.importCAD(path);
    surface.is_watertight = geo.repairManifold(surface);
    auto report = geo.analyzeManifold(surface);

    std::cout << T->label << "  Vertices          " << T->reset << surface.vertices.size() << "\n";
    std::cout << T->label << "  Faces             " << T->reset << surface.faces.size() << "\n";
    std::cout << T->label << "  Watertight        " << T->reset
              << (surface.is_watertight ? (std::string(T->ok) + "YES") : (std::string(T->err) + "NO"))
              << T->reset << "\n";
    std::cout << T->label << "  Boundary edges    " << T->reset << report.boundary_edges << "\n";
    std::cout << T->label << "  Non-manifold edges" << T->reset << " " << report.non_manifold_edges << "\n";
    std::cout << T->label << "  Isolated vertices " << T->reset << report.isolated_vertices << "\n";

    if (report.is_manifold())
        std::cout << T->ok << "  Mesh is CLEAN — ready for slicing and export.\n" << T->reset;
    else
        std::cout << T->warn << "  WARNING: issues found. Run 'heal' first.\n" << T->reset;

    return 0;
}

// ---------------------------------------------------------------------------
// heal
// ---------------------------------------------------------------------------
static int cmd_heal(const std::string& input, const std::string& output) {
    std::cout << T->label << "Healing: " << T->reset << input << "\n";
    auto surface = loadAndHeal(input);
    std::cout << T->label << "  Vertices after heal: " << T->reset << surface.vertices.size() << "\n";
    std::cout << T->label << "  Faces after heal   : " << T->reset << surface.faces.size() << "\n";
    std::cout << T->label << "  Watertight         : " << T->reset
              << (surface.is_watertight ? (std::string(T->ok) + "YES") : (std::string(T->err) + "NO"))
              << T->reset << "\n";

    cfd::makeExporter("stl_binary")->write(surface, output);
    std::cout << T->ok << "Saved: " << T->reset << output << "\n";
    return 0;
}

// ---------------------------------------------------------------------------
// slice
// ---------------------------------------------------------------------------
static int cmd_slice(const std::string& input, double z, const std::string& outPrefix) {
    std::cout << T->label << "Loading: " << T->reset << input << "\n";
    auto surface = loadAndHeal(input, true);

    std::cout << T->label << "Slicing at Z = " << T->reset << z << " ...\n";
    cfd::Slicer slicer;
    auto result = slicer.sliceAtZ(surface, z);

    if (result.contours.empty()) {
        double zmin = surface.vertices[0].pos[2], zmax = zmin;
        for (const auto& v : surface.vertices) {
            zmin = std::min(zmin, v.pos[2]);
            zmax = std::max(zmax, v.pos[2]);
        }
        std::cout << T->err << "No contours found at Z=" << z << T->reset << "\n";
        std::cout << T->dim << "  Mesh Z range: [" << zmin << ", " << zmax << "]\n" << T->reset;
        return 1;
    }

    std::size_t total = 0;
    for (const auto& c : result.contours) total += c.points.size();
    std::cout << T->label << "  Contours: " << T->reset << result.contours.size() << "\n";
    std::cout << T->label << "  Points  : " << T->reset << total << "\n";

    std::vector<cfd::SliceResult> slices{std::move(result)};
    slicer.exportSVG(slices,       outPrefix + ".svg");
    slicer.exportPolylines(slices, outPrefix + ".txt");
    std::cout << T->ok << "Saved: " << T->reset << outPrefix << ".svg\n";
    std::cout << T->ok << "Saved: " << T->reset << outPrefix << ".txt\n";
    return 0;
}

// ---------------------------------------------------------------------------
// curvature
// ---------------------------------------------------------------------------
static int cmd_curvature(const std::string& input, const std::string& output) {
    std::cout << T->label << "Loading: " << T->reset << input << "\n";
    auto surface = loadAndHeal(input, true);

    std::cout << T->label << "Computing curvature for " << T->reset
              << surface.vertices.size() << " vertices...\n";
    auto result = cfd::CurvatureAnalyzer().compute(surface);

    const auto& H = result.mean_curvature;
    double hmin = *std::min_element(H.begin(), H.end());
    double hmax = *std::max_element(H.begin(), H.end());
    std::cout << T->label << "  Min H: " << T->reset << hmin << "\n";
    std::cout << T->label << "  Max H: " << T->reset << hmax << "\n";

    cfd::makeExporter("ply", nullptr, &result)->write(surface, output);
    std::cout << T->ok << "Saved: " << T->reset << output
              << T->dim << "  (open in MeshLab or Blender to see colours)\n" << T->reset;
    return 0;
}

// ---------------------------------------------------------------------------
// export
// ---------------------------------------------------------------------------
static int cmd_export(const std::string& input, const std::string& format, const std::string& output) {
    std::cout << T->label << "Loading: " << T->reset << input << "\n";
    auto surface = loadAndHeal(input, true);

    try {
        if (format == "obj") {
            cfd::GeometryEngine geo;
            auto groups = geo.identifyFlatFaces(surface);
            cfd::makeExporter("obj", &groups)->write(surface, output);
        } else {
            cfd::makeExporter(format)->write(surface, output);
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << T->err << "Unknown format '" << format
                  << "'. Use: stl  stl_ascii  obj  ply  step\n" << T->reset;
        return 1;
    }

    std::cout << T->ok << "Saved: " << T->reset << output << "\n";
    return 0;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Theme must be set before any output
    if (hasFlag(argc, argv, "--light")) T = &LIGHT;
    else                                T = &DARK;

    if (argc < 2 || hasFlag(argc, argv, "--help") || hasFlag(argc, argv, "-h")) {
        printHelp();
        return 0;
    }

    // Skip past theme flags to find the command
    std::string command;
    std::string input;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--dark" || a == "--light") continue;
        if (command.empty()) { command = a; continue; }
        if (input.empty())   { input   = a; break; }
    }

    banner();

    if (command == "info") {
        if (input.empty()) { std::cerr << T->err << "Usage: meshmaster info <input.stl>\n" << T->reset; return 1; }
        return cmd_info(input);
    }
    if (command == "heal") {
        if (input.empty()) { std::cerr << T->err << "Usage: meshmaster heal <input.stl> [-o out.stl]\n" << T->reset; return 1; }
        return cmd_heal(input, getArg(argc, argv, "-o", "healed.stl"));
    }
    if (command == "slice") {
        if (input.empty()) { std::cerr << T->err << "Usage: meshmaster slice <input.stl> -z <value>\n" << T->reset; return 1; }
        std::string zs = getArg(argc, argv, "-z");
        if (zs.empty()) { std::cerr << T->err << "Missing -z <value>\n" << T->reset; return 1; }
        return cmd_slice(input, std::stod(zs), getArg(argc, argv, "-o", "slice"));
    }
    if (command == "curvature") {
        if (input.empty()) { std::cerr << T->err << "Usage: meshmaster curvature <input.stl> [-o out.ply]\n" << T->reset; return 1; }
        return cmd_curvature(input, getArg(argc, argv, "-o", "curvature.ply"));
    }
    if (command == "export") {
        if (input.empty()) { std::cerr << T->err << "Usage: meshmaster export <input.stl> -f <format>\n" << T->reset; return 1; }
        std::string fmt = getArg(argc, argv, "-f");
        if (fmt.empty()) { std::cerr << T->err << "Missing -f <format>\n" << T->reset; return 1; }
        std::string ext = (fmt == "stl_ascii" || fmt == "stl_binary") ? "stl" : fmt;
        return cmd_export(input, fmt, getArg(argc, argv, "-o", "surface." + ext));
    }

    std::cerr << T->err << "Unknown command '" << command << "'. Run 'meshmaster --help'\n" << T->reset;
    return 1;
}
