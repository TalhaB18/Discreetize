#pragma once
#include <string>
#include <vector>
#include <array>
#include <functional>
#include "types.hpp"

namespace cfd {

using MCV3 = std::array<double,3>;

// Run Marching Cubes over an SDF in the given axis-aligned bounding box.
// Returns a flat triangle soup: every 3 consecutive MCV3s form one triangle.
std::vector<MCV3> marching_cubes_run(
    std::function<double(double,double,double)> sdf,
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    int res);

// Convert a flat triangle soup (3 MCV3s per triangle) to a SurfaceManifold.
// Vertices are NOT merged — call GeometryEngine::mergeVertices() afterwards.
SurfaceManifold tris_to_manifold(const std::vector<MCV3>& tris);


// Box centered at (cx,cy,cz) with half-sizes (hx,hy,hz)
std::string generate_box_stl(double cx, double cy, double cz,
                              double hx, double hy, double hz);

// Ellipsoid centered at (cx,cy,cz) with radii (rx,ry,rz)
std::string generate_ellipsoid_stl(double cx, double cy, double cz,
                                   double rx, double ry, double rz,
                                   int res = 32);

// Boolean subtraction: box minus ellipsoid — returns surface mesh as binary STL
std::string boolean_subtract(double bcx, double bcy, double bcz,
                              double bhx, double bhy, double bhz,
                              double scx, double scy, double scz,
                              double srx, double sry, double srz,
                              int resolution = 80);

// Cylinder centered at (cx,cy,cz), Y-axis aligned, with given radius and height
std::string cylinder_stl(double cx, double cy, double cz,
                          double radius, double height, int segs = 32);

// Cone with apex at top (cx, cy+height/2, cz) and base at bottom
std::string cone_stl(double cx, double cy, double cz,
                     double radius, double height, int segs = 24);

// Torus centered at (cx,cy,cz) around Y axis
std::string torus_stl(double cx, double cy, double cz,
                       double major_r, double minor_r,
                       int major_segs = 32, int minor_segs = 16);

// Flat horizontal (XZ) plane subdivided into subdiv x subdiv quads
std::string plane_stl(double cx, double cy, double cz,
                       double width, double depth, int subdiv = 4);

// Icosphere: subdivided icosahedron, uniform triangles, subdivisions 0-4
std::string icosphere_stl(double cx, double cy, double cz,
                           double radius, int subdivisions = 2);

} // namespace cfd
