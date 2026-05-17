#include "mesh_sdf.hpp"
#include <cmath>
#include <limits>
#include <array>

namespace cfd {

// ---------------------------------------------------------------------------
// Vector helpers (local, not exported)
// ---------------------------------------------------------------------------
static inline double v3dot(const Vec3d& a, const Vec3d& b){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static inline Vec3d v3sub(const Vec3d& a, const Vec3d& b){
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
static inline Vec3d v3cross(const Vec3d& a, const Vec3d& b){
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
static inline double v3norm2(const Vec3d& a){ return v3dot(a,a); }
static inline double v3norm(const Vec3d& a){ return std::sqrt(v3norm2(a)); }

// Closest-point-on-triangle distance² using Eberly's method.
static double ptTriDistSq(const Vec3d& p,
                           const Vec3d& A, const Vec3d& B, const Vec3d& C)
{
    Vec3d AB = v3sub(B,A), AC = v3sub(C,A), AP = v3sub(p,A);
    double d1=v3dot(AB,AP), d2=v3dot(AC,AP);
    if(d1<=0&&d2<=0) return v3norm2(AP);

    Vec3d BP = v3sub(p,B);
    double d3=v3dot(AB,BP), d4=v3dot(AC,BP);
    if(d3>=0&&d4<=d3){ return v3norm2(BP); }

    double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0){
        double v=d1/(d1-d3);
        Vec3d q={A[0]+v*AB[0],A[1]+v*AB[1],A[2]+v*AB[2]};
        return v3norm2(v3sub(p,q));
    }

    Vec3d CP = v3sub(p,C);
    double d5=v3dot(AB,CP), d6=v3dot(AC,CP);
    if(d6>=0&&d5<=d6){ return v3norm2(CP); }

    double vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0){
        double w=d2/(d2-d6);
        Vec3d q={A[0]+w*AC[0],A[1]+w*AC[1],A[2]+w*AC[2]};
        return v3norm2(v3sub(p,q));
    }

    double va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){
        double w=(d4-d3)/((d4-d3)+(d5-d6));
        Vec3d BC = v3sub(C,B);
        Vec3d q={B[0]+w*BC[0],B[1]+w*BC[1],B[2]+w*BC[2]};
        return v3norm2(v3sub(p,q));
    }

    double denom=1.0/(va+vb+vc);
    double v=vb*denom, w=vc*denom;
    Vec3d q={A[0]+v*AB[0]+w*AC[0],A[1]+v*AB[1]+w*AC[1],A[2]+v*AB[2]+w*AC[2]};
    return v3norm2(v3sub(p,q));
}

// AABB min squared distance to point.
static double aabbMinDistSq(const AABB& box, const Vec3d& p){
    double d=0;
    for(int k=0;k<3;k++){
        if(p[k]<box.lo[k]){ double t=box.lo[k]-p[k]; d+=t*t; }
        else if(p[k]>box.hi[k]){ double t=p[k]-box.hi[k]; d+=t*t; }
    }
    return d;
}

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
MeshSDF::MeshSDF(const SurfaceManifold& mesh) : mesh_(mesh) {
    bvh_.build(mesh);

    const double kInf = std::numeric_limits<double>::infinity();
    bbox_lo_ = { kInf, kInf, kInf};
    bbox_hi_ = {-kInf,-kInf,-kInf};
    for(const auto& v : mesh.vertices){
        for(int k=0;k<3;k++){
            if(v.pos[k]<bbox_lo_[k]) bbox_lo_[k]=v.pos[k];
            if(v.pos[k]>bbox_hi_[k]) bbox_hi_[k]=v.pos[k];
        }
    }
    // small padding
    double diag=v3norm(v3sub(bbox_hi_,bbox_lo_));
    for(int k=0;k<3;k++){
        bbox_lo_[k] -= 0.01*diag;
        bbox_hi_[k] += 0.01*diag;
    }
}

// ---------------------------------------------------------------------------
// eval — unsigned distance then sign from ray casting
// ---------------------------------------------------------------------------
double MeshSDF::eval(const Vec3d& p) const {
    if(mesh_.faces.empty()) return 1.0;
    double dist=0;
    uint32_t tri_idx = bvh_.nearestTriangle(p, dist);
    (void)tri_idx;
    return isInside(p) ? -dist : dist;
}

// ---------------------------------------------------------------------------
// isInside — majority vote over 3 axis-aligned ray casts
// ---------------------------------------------------------------------------
bool MeshSDF::isInside(const Vec3d& p) const {
    int votes = 0;
    Vec3d dirs[3] = {{1,0,0},{0,1,0},{0,0,1}};
    // Offset slightly to avoid hitting edges/vertices exactly
    Vec3d po = {p[0]+1.3e-7, p[1]+1.7e-7, p[2]+2.1e-7};
    for(const auto& d : dirs){
        if((raycastCount(po, d) % 2) == 1) votes++;
    }
    return votes >= 2;
}

// ---------------------------------------------------------------------------
// raycastCount — count ray-triangle intersections (Möller–Trumbore)
// ---------------------------------------------------------------------------
int MeshSDF::raycastCount(const Vec3d& orig, const Vec3d& dir) const {
    constexpr double EPS = 1e-10;
    int count = 0;
    for(const auto& tri : mesh_.faces){
        const Vec3d& v0 = mesh_.vertices[tri.vi[0]].pos;
        const Vec3d& v1 = mesh_.vertices[tri.vi[1]].pos;
        const Vec3d& v2 = mesh_.vertices[tri.vi[2]].pos;
        Vec3d e1=v3sub(v1,v0), e2=v3sub(v2,v0);
        Vec3d h=v3cross(dir,e2);
        double a=v3dot(e1,h);
        if(std::abs(a)<EPS) continue;
        double f=1.0/a;
        Vec3d s=v3sub(orig,v0);
        double u=f*v3dot(s,h);
        if(u<0||u>1) continue;
        Vec3d q=v3cross(s,e1);
        double v=f*v3dot(dir,q);
        if(v<0||u+v>1) continue;
        double t=f*v3dot(e2,q);
        if(t>EPS) count++;
    }
    return count;
}

} // namespace cfd
