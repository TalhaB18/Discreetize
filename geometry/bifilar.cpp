#include "bifilar.hpp"
#include "primitives.hpp"
#include "geometry_engine.hpp"
#include <cmath>
#include <algorithm>

namespace cfd {

std::function<double(double,double,double)> makeBifilarSDF(const BifilarConfig& cfg)
{
    return [cfg](double x, double y, double z) -> double {
        const double omega = (2.0 * M_PI * cfg.revolutions) / cfg.H;
        const double theta = omega * z;
        const double cx1 =  cfg.offset * std::cos(theta);
        const double cy1 =  cfg.offset * std::sin(theta);
        const double cx2 = -cfg.offset * std::cos(theta);
        const double cy2 = -cfg.offset * std::sin(theta);

        // Outer cylinder: negative inside
        double sdf_outer = std::sqrt(x*x + y*y) - cfg.R;

        // Inner cylinders: negative inside each
        double sdf_i1 = std::sqrt((x-cx1)*(x-cx1)+(y-cy1)*(y-cy1)) - cfg.r;
        double sdf_i2 = std::sqrt((x-cx2)*(x-cx2)+(y-cy2)*(y-cy2)) - cfg.r;

        // Z caps
        double sdf_zlo = -z;        // positive below z=0
        double sdf_zhi = z - cfg.H; // positive above z=H

        // Fluid domain: inside outer AND outside both inners AND 0<=z<=H
        // SDF < 0 inside fluid:
        return std::max({sdf_outer, -sdf_i1, -sdf_i2, sdf_zlo, sdf_zhi});
    };
}

SurfaceManifold generateBifilarSurface(const BifilarConfig& cfg)
{
    auto sdf = makeBifilarSDF(cfg);

    const double pad = 0.02 * cfg.R;
    double x0 = -(cfg.R + pad), x1 = cfg.R + pad;
    double y0 = -(cfg.R + pad), y1 = cfg.R + pad;
    double z0 = -pad,            z1 = cfg.H + pad;

    // Isotropic voxel resolution
    int res = std::max(20, static_cast<int>(cfg.resolution * 2.0 * (cfg.R + pad)));

    auto tris = marching_cubes_run(sdf, x0, y0, z0, x1, y1, z1, res);

    auto surface = tris_to_manifold(tris);

    GeometryEngine geo;
    geo.mergeVertices(surface, 1e-8);
    (void)geo.repairManifold(surface);

    // Tag all faces as Wall — required by MeshEngine::inflatePrisms()
    BoundaryPatch wall;
    wall.name = "wall";
    wall.type = BoundaryType::Wall;
    wall.face_ids.resize(surface.faces.size());
    for(uint32_t i = 0; i < (uint32_t)surface.faces.size(); ++i)
        wall.face_ids[i] = i;
    surface.patches.push_back(std::move(wall));

    return surface;
}

SurfaceManifold generateBifilarInnerCylinder(const BifilarConfig& cfg, int which)
{
    const double omega = (2.0 * M_PI * cfg.revolutions) / cfg.H;
    const double sign  = (which == 1) ? 1.0 : -1.0;

    auto sdf = [cfg, omega, sign](double x, double y, double z) -> double {
        const double theta = omega * z;
        const double cx = sign * cfg.offset * std::cos(theta);
        const double cy = sign * cfg.offset * std::sin(theta);
        const double dx = x - cx, dy = y - cy;
        return std::max({ std::sqrt(dx*dx + dy*dy) - cfg.r, -z, z - cfg.H });
    };

    const double pad = 0.02 * cfg.R;
    double x0 = -(cfg.R + pad), x1 = cfg.R + pad;
    double y0 = -(cfg.R + pad), y1 = cfg.R + pad;
    double z0 = -pad,           z1 = cfg.H + pad;
    int res = std::max(20, static_cast<int>(cfg.resolution * 2.0 * (cfg.R + pad)));

    auto tris    = marching_cubes_run(sdf, x0, y0, z0, x1, y1, z1, res);
    auto surface = tris_to_manifold(tris);

    GeometryEngine geo;
    geo.mergeVertices(surface, 1e-8);
    (void)geo.repairManifold(surface);

    BoundaryPatch wall;
    wall.name = "wall"; wall.type = BoundaryType::Wall;
    wall.face_ids.resize(surface.faces.size());
    for (uint32_t i = 0; i < (uint32_t)surface.faces.size(); ++i) wall.face_ids[i] = i;
    surface.patches.push_back(std::move(wall));

    return surface;
}

} // namespace cfd
