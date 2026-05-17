#include "primitives.hpp"
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <map>
#include <vector>
#include <algorithm>

namespace cfd {

using V3 = std::array<double,3>;

// ─── Binary STL helpers ───────────────────────────────────────────────────────

static void wf  (std::string& s, float    v){ char b[4]; std::memcpy(b,&v,4); s.append(b,4); }
static void wu32(std::string& s, uint32_t v){ char b[4]; std::memcpy(b,&v,4); s.append(b,4); }
static void wu16(std::string& s, uint16_t v){ char b[2]; std::memcpy(b,&v,2); s.append(b,2); }

static V3 vsub  (const V3& a,const V3& b){ return{a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
static V3 vcross(const V3& a,const V3& b){ return{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}; }
static double vdot(const V3& a,const V3& b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static V3 vnorm(const V3& a){
    double l=std::sqrt(vdot(a,a));
    return l<1e-12?V3{0,0,1}:V3{a[0]/l,a[1]/l,a[2]/l};
}

static void write_tri(std::string& out,const V3& a,const V3& b,const V3& c){
    V3 n=vnorm(vcross(vsub(b,a),vsub(c,a)));
    wf(out,(float)n[0]); wf(out,(float)n[1]); wf(out,(float)n[2]);
    for(auto& v:{a,b,c}){ wf(out,(float)v[0]); wf(out,(float)v[1]); wf(out,(float)v[2]); }
    wu16(out,0);
}

static std::string make_stl(const std::vector<V3>& verts){
    uint32_t n=(uint32_t)(verts.size()/3);
    std::string out; out.reserve(84+n*50);
    out.append(80,'\0');
    wu32(out,n);
    for(uint32_t i=0;i<n;++i) write_tri(out,verts[i*3],verts[i*3+1],verts[i*3+2]);
    return out;
}

// Feature-preserving Laplacian smoothing on a triangle soup (in-place).
// Vertices on edges where adjacent face normals differ by >= sharp_deg are
// locked and never moved, so box corners/edges stay crisp while curved
// surfaces (sphere crater) become smooth.
static void smooth_mesh(std::vector<V3>& tris, int iters = 4, double alpha = 0.5,
                        double sharp_deg = 30.0)
{
    if (tris.empty()) return;

    // ── 1. Build unique vertex index (exact float equality — MC verts share coords) ──
    std::map<std::tuple<float,float,float>, int> vidx;
    std::vector<V3>  uverts;
    std::vector<int> triIdx(tris.size());
    for (int i = 0; i < (int)tris.size(); i++) {
        auto key = std::make_tuple((float)tris[i][0],(float)tris[i][1],(float)tris[i][2]);
        auto [it, ins] = vidx.emplace(key, (int)uverts.size());
        if (ins) uverts.push_back(tris[i]);
        triIdx[i] = it->second;
    }
    int nV   = (int)uverts.size();
    int nTri = (int)tris.size() / 3;

    // ── 2. Per-triangle face normals ──────────────────────────────────────────
    std::vector<V3> faceN(nTri);
    for (int t = 0; t < nTri; t++) {
        V3 A=uverts[triIdx[t*3]], B=uverts[triIdx[t*3+1]], C=uverts[triIdx[t*3+2]];
        faceN[t] = vnorm(vcross(vsub(B,A), vsub(C,A)));
    }

    // ── 3. Edge → triangle list (directed edge stored as min/max pair) ────────
    std::map<std::pair<int,int>, std::vector<int>> edgeTris;
    for (int t = 0; t < nTri; t++) {
        int a=triIdx[t*3], b=triIdx[t*3+1], c=triIdx[t*3+2];
        for (auto [u,v] : std::initializer_list<std::pair<int,int>>{{a,b},{b,c},{a,c}}) {
            if (u > v) std::swap(u, v);
            edgeTris[{u,v}].push_back(t);
        }
    }

    // ── 4. Mark sharp vertices ────────────────────────────────────────────────
    const double cosT = std::cos(sharp_deg * M_PI / 180.0);
    std::vector<bool> sharp(nV, false);
    for (auto& [edge, ts] : edgeTris) {
        for (int i = 0; i < (int)ts.size(); i++) {
            for (int j = i+1; j < (int)ts.size(); j++) {
                if (vdot(faceN[ts[i]], faceN[ts[j]]) < cosT) {
                    sharp[edge.first]  = true;
                    sharp[edge.second] = true;
                }
            }
        }
    }

    // ── 5. Build neighbor lists ───────────────────────────────────────────────
    std::vector<std::vector<int>> nbrs(nV);
    for (int t = 0; t < nTri; t++) {
        int a=triIdx[t*3], b=triIdx[t*3+1], c=triIdx[t*3+2];
        nbrs[a].push_back(b); nbrs[a].push_back(c);
        nbrs[b].push_back(a); nbrs[b].push_back(c);
        nbrs[c].push_back(a); nbrs[c].push_back(b);
    }

    // ── 6. Laplacian iterations — skip locked (sharp) vertices ───────────────
    std::vector<V3> nv = uverts;
    for (int iter = 0; iter < iters; iter++) {
        for (int i = 0; i < nV; i++) {
            if (sharp[i] || nbrs[i].empty()) continue;
            double sx=0, sy=0, sz=0;
            double n = (double)nbrs[i].size();
            for (int j : nbrs[i]) { sx+=uverts[j][0]; sy+=uverts[j][1]; sz+=uverts[j][2]; }
            nv[i] = { uverts[i][0]*(1-alpha) + sx/n*alpha,
                      uverts[i][1]*(1-alpha) + sy/n*alpha,
                      uverts[i][2]*(1-alpha) + sz/n*alpha };
        }
        uverts = nv;
    }

    // ── 7. Write back to flat triangle soup ───────────────────────────────────
    for (int i = 0; i < (int)tris.size(); i++) tris[i] = uverts[triIdx[i]];
}

// ─── Box (analytical) ────────────────────────────────────────────────────────

std::string generate_box_stl(double cx,double cy,double cz,
                              double hx,double hy,double hz){
    std::vector<V3> t; t.reserve(36);
    auto q=[&](V3 a,V3 b,V3 c,V3 d){
        t.push_back(a);t.push_back(b);t.push_back(c);
        t.push_back(a);t.push_back(c);t.push_back(d);
    };
    double x0=cx-hx,x1=cx+hx,y0=cy-hy,y1=cy+hy,z0=cz-hz,z1=cz+hz;
    q({x0,y0,z0},{x0,y0,z1},{x0,y1,z1},{x0,y1,z0}); // -X
    q({x1,y0,z0},{x1,y1,z0},{x1,y1,z1},{x1,y0,z1}); // +X
    q({x0,y0,z0},{x1,y0,z0},{x1,y0,z1},{x0,y0,z1}); // -Y
    q({x0,y1,z0},{x0,y1,z1},{x1,y1,z1},{x1,y1,z0}); // +Y
    q({x0,y0,z0},{x0,y1,z0},{x1,y1,z0},{x1,y0,z0}); // -Z
    q({x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1}); // +Z
    return make_stl(t);
}

// ─── Ellipsoid (UV) ───────────────────────────────────────────────────────────

std::string generate_ellipsoid_stl(double cx,double cy,double cz,
                                   double rx,double ry,double rz,int res){
    const double pi=std::acos(-1.0);
    std::vector<V3> t;
    for(int lat=0;lat<res;++lat){
        double t0=pi*lat/res, t1=pi*(lat+1)/res;
        for(int lon=0;lon<res*2;++lon){
            double p0=2*pi*lon/(res*2), p1=2*pi*(lon+1)/(res*2);
            auto v=[&](double th,double ph)->V3{
                return{cx+rx*std::sin(th)*std::cos(ph),
                       cy+ry*std::cos(th),
                       cz+rz*std::sin(th)*std::sin(ph)};
            };
            V3 v00=v(t0,p0),v01=v(t0,p1),v10=v(t1,p0),v11=v(t1,p1);
            if(lat==0){
                // North-pole cap: single triangle per longitude slice
                t.push_back(v00);t.push_back(v11);t.push_back(v10);
            } else if(lat==res-1){
                // South-pole cap: single triangle per longitude slice
                t.push_back(v00);t.push_back(v01);t.push_back(v10);
            } else {
                // Middle band: quad as two triangles
                t.push_back(v00);t.push_back(v10);t.push_back(v11);
                t.push_back(v00);t.push_back(v11);t.push_back(v01);
            }
        }
    }
    return make_stl(t);
}

// ─── Marching Cubes tables ────────────────────────────────────────────────────

static const int MC_EDGE[256]={
    0x000,0x109,0x203,0x30a,0x406,0x50f,0x605,0x70c,
    0x80c,0x905,0xa0f,0xb06,0xc0a,0xd03,0xe09,0xf00,
    0x190,0x099,0x393,0x29a,0x596,0x49f,0x795,0x69c,
    0x99c,0x895,0xb9f,0xa96,0xd9a,0xc93,0xf99,0xe90,
    0x230,0x339,0x033,0x13a,0x636,0x73f,0x435,0x53c,
    0xa3c,0xb35,0x83f,0x936,0xe3a,0xf33,0xc39,0xd30,
    0x3a0,0x2a9,0x1a3,0x0aa,0x7a6,0x6af,0x5a5,0x4ac,
    0xbac,0xaa5,0x9af,0x8a6,0xfaa,0xea3,0xda9,0xca0,
    0x460,0x569,0x663,0x76a,0x066,0x16f,0x265,0x36c,
    0xc6c,0xd65,0xe6f,0xf66,0x86a,0x963,0xa69,0xb60,
    0x5f0,0x4f9,0x7f3,0x6fa,0x1f6,0x0ff,0x3f5,0x2fc,
    0xdfc,0xcf5,0xfff,0xef6,0x9fa,0x8f3,0xbf9,0xaf0,
    0x650,0x759,0x453,0x55a,0x256,0x35f,0x055,0x15c,
    0xe5c,0xf55,0xc5f,0xd56,0xa5a,0xb53,0x859,0x950,
    0x7c0,0x6c9,0x5c3,0x4ca,0x3c6,0x2cf,0x1c5,0x0cc,
    0xfcc,0xec5,0xdcf,0xcc6,0xbca,0xac3,0x9c9,0x8c0,
    0x8c0,0x9c9,0xac3,0xbca,0xcc6,0xdcf,0xec5,0xfcc,
    0x0cc,0x1c5,0x2cf,0x3c6,0x4ca,0x5c3,0x6c9,0x7c0,
    0x950,0x859,0xb53,0xa5a,0xd56,0xc5f,0xf55,0xe5c,
    0x15c,0x055,0x35f,0x256,0x55a,0x453,0x759,0x650,
    0xaf0,0xbf9,0x8f3,0x9fa,0xef6,0xfff,0xcf5,0xdfc,
    0x2fc,0x3f5,0x0ff,0x1f6,0x6fa,0x7f3,0x4f9,0x5f0,
    0xb60,0xa69,0x963,0x86a,0xf66,0xe6f,0xd65,0xc6c,
    0x36c,0x265,0x16f,0x066,0x76a,0x663,0x569,0x460,
    0xca0,0xda9,0xea3,0xfaa,0x8a6,0x9af,0xaa5,0xbac,
    0x4ac,0x5a5,0x6af,0x7a6,0x0aa,0x1a3,0x2a9,0x3a0,
    0xd30,0xc39,0xf33,0xe3a,0x936,0x83f,0xb35,0xa3c,
    0x53c,0x435,0x73f,0x636,0x13a,0x033,0x339,0x230,
    0xe90,0xf99,0xc93,0xd9a,0xa96,0xb9f,0x895,0x99c,
    0x69c,0x795,0x49f,0x596,0x29a,0x393,0x099,0x190,
    0xf00,0xe09,0xd03,0xc0a,0xb06,0xa0f,0x905,0x80c,
    0x70c,0x605,0x50f,0x406,0x30a,0x203,0x109,0x000
};

static const int MC_TRIS[256][16]={
    {-1},{0,8,3,-1},{0,1,9,-1},{1,8,3,9,8,1,-1},
    {1,2,10,-1},{0,8,3,1,2,10,-1},{9,2,10,0,2,9,-1},{2,8,3,2,10,8,10,9,8,-1},
    {3,11,2,-1},{0,11,2,8,11,0,-1},{1,9,0,2,3,11,-1},{1,11,2,1,9,11,9,8,11,-1},
    {3,10,1,11,10,3,-1},{0,10,1,0,8,10,8,11,10,-1},{3,9,0,3,11,9,11,10,9,-1},{9,8,10,10,8,11,-1},
    {4,7,8,-1},{4,3,0,7,3,4,-1},{0,1,9,8,4,7,-1},{4,1,9,4,7,1,7,3,1,-1},
    {1,2,10,8,4,7,-1},{3,4,7,3,0,4,1,2,10,-1},{9,2,10,9,0,2,8,4,7,-1},{2,10,9,2,9,7,2,7,3,7,9,4,-1},
    {8,4,7,3,11,2,-1},{11,4,7,11,2,4,2,0,4,-1},{9,0,1,8,4,7,2,3,11,-1},{4,7,11,9,4,11,9,11,2,9,2,1,-1},
    {3,10,1,3,11,10,7,8,4,-1},{1,11,10,1,4,11,1,0,4,7,11,4,-1},{4,7,8,9,0,11,9,11,10,11,0,3,-1},{4,7,11,4,11,9,9,11,10,-1},
    {9,5,4,-1},{9,5,4,0,8,3,-1},{0,5,4,1,5,0,-1},{8,5,4,8,3,5,3,1,5,-1},
    {1,2,10,9,5,4,-1},{3,0,8,1,2,10,4,9,5,-1},{5,2,10,5,4,2,4,0,2,-1},{2,10,5,3,2,5,3,5,4,3,4,8,-1},
    {9,5,4,2,3,11,-1},{0,11,2,0,8,11,4,9,5,-1},{0,5,4,0,1,5,2,3,11,-1},{2,1,5,2,5,8,2,8,11,4,8,5,-1},
    {10,3,11,10,1,3,9,5,4,-1},{4,9,5,0,8,1,8,10,1,8,11,10,-1},{5,4,0,5,0,11,5,11,10,11,0,3,-1},{5,4,8,5,8,10,10,8,11,-1},
    {9,7,8,5,7,9,-1},{9,3,0,9,5,3,5,7,3,-1},{0,7,8,0,1,7,1,5,7,-1},{1,5,3,3,5,7,-1},
    {9,7,8,9,5,7,10,1,2,-1},{10,1,2,9,5,0,5,3,0,5,7,3,-1},{8,0,2,8,2,5,8,5,7,10,5,2,-1},{2,10,5,2,5,3,3,5,7,-1},
    {7,9,5,7,8,9,3,11,2,-1},{9,5,7,9,7,2,9,2,0,2,7,11,-1},{2,3,11,0,1,8,1,7,8,1,5,7,-1},{11,2,1,11,1,7,7,1,5,-1},
    {9,5,8,8,5,7,10,1,3,10,3,11,-1},{5,7,0,5,0,9,7,11,0,1,0,10,11,10,0,-1},{11,10,0,11,0,3,10,5,0,8,0,7,5,7,0,-1},{11,10,5,7,11,5,-1},
    {10,6,5,-1},{0,8,3,5,10,6,-1},{9,0,1,5,10,6,-1},{1,8,3,1,9,8,5,10,6,-1},
    {1,6,5,2,6,1,-1},{1,6,5,1,2,6,3,0,8,-1},{9,6,5,9,0,6,0,2,6,-1},{5,9,8,5,8,2,5,2,6,3,2,8,-1},
    {2,3,11,10,6,5,-1},{11,0,8,11,2,0,10,6,5,-1},{0,1,9,2,3,11,5,10,6,-1},{5,10,6,1,9,2,9,11,2,9,8,11,-1},
    {6,3,11,6,5,3,5,1,3,-1},{0,8,11,0,11,5,0,5,1,5,11,6,-1},{3,11,6,0,3,6,0,6,5,0,5,9,-1},{6,5,9,6,9,11,11,9,8,-1},
    {5,10,6,4,7,8,-1},{4,3,0,4,7,3,6,5,10,-1},{1,9,0,5,10,6,8,4,7,-1},{10,6,5,1,9,7,1,7,3,7,9,4,-1},
    {6,1,2,6,5,1,4,7,8,-1},{1,2,5,5,2,6,3,0,4,3,4,7,-1},{8,4,7,9,0,5,0,6,5,0,2,6,-1},{7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,-1},
    {3,11,2,7,8,4,10,6,5,-1},{5,10,6,4,7,2,4,2,0,2,7,11,-1},{0,1,9,4,7,8,2,3,11,5,10,6,-1},{9,2,1,9,11,2,9,4,11,7,11,4,5,10,6,-1},
    {8,4,7,3,11,5,3,5,1,5,11,6,-1},{5,1,11,5,11,6,1,0,11,7,11,4,0,4,11,-1},{0,5,9,0,6,5,0,3,6,11,6,3,8,4,7,-1},{6,5,9,6,9,11,4,7,9,7,11,9,-1},
    {10,4,9,6,4,10,-1},{4,10,6,4,9,10,0,8,3,-1},{10,0,1,10,6,0,6,4,0,-1},{8,3,1,8,1,6,8,6,4,6,1,10,-1},
    {1,4,9,1,2,4,2,6,4,-1},{3,0,8,1,2,9,2,4,9,2,6,4,-1},{0,2,4,4,2,6,-1},{8,3,2,8,2,4,4,2,6,-1},
    {10,4,9,10,6,4,11,2,3,-1},{0,8,2,2,8,11,4,9,10,4,10,6,-1},{3,11,2,0,1,6,0,6,4,6,1,10,-1},{6,4,1,6,1,10,4,8,1,2,1,11,8,11,1,-1},
    {9,6,4,9,3,6,9,1,3,11,6,3,-1},{8,11,1,8,1,0,11,6,1,9,1,4,6,4,1,-1},{3,11,6,3,6,0,0,6,4,-1},{6,4,8,11,6,8,-1},
    {7,10,6,7,8,10,8,9,10,-1},{0,7,3,0,10,7,0,9,10,6,7,10,-1},{10,6,7,1,10,7,1,7,8,1,8,0,-1},{10,6,7,10,7,1,1,7,3,-1},
    {1,2,6,1,6,8,1,8,9,8,6,7,-1},{2,6,9,2,9,1,6,7,9,0,9,3,7,3,9,-1},{7,8,0,7,0,6,6,0,2,-1},{7,3,2,6,7,2,-1},
    {2,3,11,10,6,8,10,8,9,8,6,7,-1},{2,0,7,2,7,11,0,9,7,6,7,10,9,10,7,-1},{1,8,0,1,7,8,1,10,7,6,7,10,2,3,11,-1},{11,2,1,11,1,7,10,6,1,6,7,1,-1},
    {8,9,1,8,1,3,9,6,1,11,1,7,6,7,1,-1},{0,9,1,11,6,7,-1},{7,8,0,7,0,6,3,11,0,11,6,0,-1},{7,11,6,-1},
    {7,6,11,-1},{3,0,8,11,7,6,-1},{0,1,9,11,7,6,-1},{8,1,9,8,3,1,11,7,6,-1},
    {10,1,2,6,11,7,-1},{1,2,10,3,0,8,6,11,7,-1},{2,9,0,2,10,9,6,11,7,-1},{6,11,7,2,10,3,10,8,3,10,9,8,-1},
    {7,2,3,6,2,7,-1},{7,0,8,7,6,0,6,2,0,-1},{2,7,6,2,3,7,0,1,9,-1},{1,6,2,1,8,6,1,9,8,8,7,6,-1},
    {10,7,6,10,1,7,1,3,7,-1},{10,7,6,1,7,10,1,8,7,1,0,8,-1},{0,3,7,0,7,10,0,10,9,6,10,7,-1},{7,6,10,7,10,8,8,10,9,-1},
    {6,8,4,11,8,6,-1},{3,6,11,3,0,6,0,4,6,-1},{8,6,11,8,4,6,9,0,1,-1},{9,4,6,9,6,3,9,3,1,11,3,6,-1},
    {6,8,4,6,11,8,2,10,1,-1},{1,2,10,3,0,11,0,6,11,0,4,6,-1},{4,11,8,4,6,11,0,2,9,2,10,9,-1},{10,9,3,10,3,2,9,4,3,11,3,6,4,6,3,-1},
    {8,2,3,8,4,2,4,6,2,-1},{0,4,2,4,6,2,-1},{1,9,0,2,3,4,2,4,6,4,3,8,-1},{1,9,4,1,4,2,2,4,6,-1},
    {8,1,3,8,6,1,8,4,6,6,10,1,-1},{10,1,0,10,0,6,6,0,4,-1},{4,6,3,4,3,8,6,10,3,0,3,9,10,9,3,-1},{10,9,4,6,10,4,-1},
    {4,9,5,7,6,11,-1},{0,8,3,4,9,5,11,7,6,-1},{5,0,1,5,4,0,7,6,11,-1},{11,7,6,8,3,4,3,5,4,3,1,5,-1},
    {9,5,4,10,1,2,7,6,11,-1},{6,11,7,1,2,10,0,8,3,4,9,5,-1},{7,6,11,5,4,10,4,2,10,4,0,2,-1},{3,4,8,3,5,4,3,2,5,10,5,2,11,7,6,-1},
    {7,2,3,7,6,2,5,4,9,-1},{9,5,4,0,8,6,0,6,2,6,8,7,-1},{3,6,2,3,7,6,1,5,0,5,4,0,-1},{6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,-1},
    {9,5,4,10,1,6,1,7,6,1,3,7,-1},{1,6,10,1,7,6,1,0,7,8,7,0,9,5,4,-1},{4,0,10,4,10,5,0,3,10,6,10,7,3,7,10,-1},{7,6,10,7,10,8,5,4,10,4,8,10,-1},
    {6,9,5,6,11,9,11,8,9,-1},{3,6,11,0,6,3,0,5,6,0,9,5,-1},{0,11,8,0,5,11,0,1,5,5,6,11,-1},{6,11,3,6,3,5,5,3,1,-1},
    {1,2,10,9,5,11,9,11,8,11,5,6,-1},{0,11,3,0,6,11,0,9,6,5,6,9,1,2,10,-1},{11,8,5,11,5,6,8,0,5,10,5,2,0,2,5,-1},{6,11,3,6,3,5,2,10,3,10,5,3,-1},
    {5,8,9,5,2,8,5,6,2,3,8,2,-1},{9,5,6,9,6,0,0,6,2,-1},{1,5,8,1,8,0,5,6,8,3,8,2,6,2,8,-1},{1,5,6,2,1,6,-1},
    {1,3,6,1,6,10,3,8,6,5,6,9,8,9,6,-1},{10,1,0,10,0,6,9,5,0,5,6,0,-1},{0,3,8,5,6,10,-1},{10,5,6,-1},
    {11,5,10,7,5,11,-1},{11,5,10,11,7,5,8,3,0,-1},{5,11,7,5,10,11,1,9,0,-1},{10,7,5,10,11,7,9,8,1,8,3,1,-1},
    {11,1,2,11,7,1,7,5,1,-1},{0,8,3,1,2,7,1,7,5,7,2,11,-1},{9,7,5,9,2,7,9,0,2,2,11,7,-1},{7,5,2,7,2,11,5,9,2,3,2,8,9,8,2,-1},
    {2,5,10,2,3,5,3,7,5,-1},{8,2,0,8,5,2,8,7,5,10,2,5,-1},{9,0,1,5,10,3,5,3,7,3,10,2,-1},{9,8,2,9,2,1,8,7,2,10,2,5,7,5,2,-1},
    {1,3,5,3,7,5,-1},{0,8,7,0,7,1,1,7,5,-1},{9,0,3,9,3,5,5,3,7,-1},{9,8,7,5,9,7,-1},
    {5,8,4,5,10,8,10,11,8,-1},{5,0,4,5,11,0,5,10,11,11,3,0,-1},{0,1,9,8,4,10,8,10,11,10,4,5,-1},{10,11,4,10,4,5,11,3,4,9,4,1,3,1,4,-1},
    {2,5,1,2,8,5,2,11,8,4,5,8,-1},{0,4,11,0,11,3,4,5,11,2,11,1,5,1,11,-1},{0,2,5,0,5,9,2,11,5,4,5,8,11,8,5,-1},{9,4,5,2,11,3,-1},
    {2,5,10,3,5,2,3,4,5,3,8,4,-1},{5,10,2,5,2,4,4,2,0,-1},{3,10,2,3,5,10,3,8,5,4,5,8,0,1,9,-1},{5,10,2,5,2,4,1,9,2,9,4,2,-1},
    {8,4,5,8,5,3,3,5,1,-1},{0,4,5,1,0,5,-1},{8,4,5,8,5,3,9,0,5,0,3,5,-1},{9,4,5,-1},
    {4,11,7,4,9,11,9,10,11,-1},{0,8,3,4,9,7,9,11,7,9,10,11,-1},{1,10,11,1,11,4,1,4,0,7,4,11,-1},{3,1,4,3,4,8,1,10,4,7,4,11,10,11,4,-1},
    {4,11,7,9,11,4,9,2,11,9,1,2,-1},{9,7,4,9,11,7,9,1,11,2,11,1,0,8,3,-1},{11,7,4,11,4,2,2,4,0,-1},{11,7,4,11,4,2,8,3,4,3,2,4,-1},
    {2,9,10,2,7,9,2,3,7,7,4,9,-1},{9,10,7,9,7,4,10,2,7,8,7,0,2,0,7,-1},{3,7,10,3,10,2,7,4,10,1,10,0,4,0,10,-1},{1,10,2,8,7,4,-1},
    {4,9,1,4,1,7,7,1,3,-1},{4,9,1,4,1,7,0,8,1,8,7,1,-1},{4,0,3,7,4,3,-1},{4,8,7,-1},
    {9,10,8,10,11,8,-1},{3,0,9,3,9,11,11,9,10,-1},{0,1,10,0,10,8,8,10,11,-1},{3,1,10,11,3,10,-1},
    {1,2,11,1,11,9,9,11,8,-1},{3,0,9,3,9,11,1,2,9,2,11,9,-1},{0,2,11,8,0,11,-1},{3,2,11,-1},
    {2,3,8,2,8,10,10,8,9,-1},{9,10,2,0,9,2,-1},{2,3,8,2,8,10,0,1,8,1,10,8,-1},{1,10,2,-1},
    {1,3,8,9,1,8,-1},{0,9,1,-1},{0,3,8,-1},{-1}
};

static const int MC_EV[12][2]={{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
static const double MC_OFF[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};

static V3 edge_interp(const V3 vp[8],const double vv[8],int e){
    int a=MC_EV[e][0],b=MC_EV[e][1];
    double denom=vv[a]-vv[b];
    double t=std::abs(denom)<1e-12 ? 0.5 : vv[a]/denom;
    t=t<0.0?0.0:t>1.0?1.0:t;
    return{vp[a][0]+t*(vp[b][0]-vp[a][0]),vp[a][1]+t*(vp[b][1]-vp[a][1]),vp[a][2]+t*(vp[b][2]-vp[a][2])};
}

static std::vector<V3> marching_cubes_verts(
    std::function<double(double,double,double)> sdf,
    double x0,double y0,double z0,double x1,double y1,double z1,int res)
{
    double dx=(x1-x0)/res,dy=(y1-y0)/res,dz=(z1-z0)/res;
    std::vector<double> grid((res+1)*(res+1)*(res+1));
    auto gidx=[&](int ix,int iy,int iz){ return iz*(res+1)*(res+1)+iy*(res+1)+ix; };
    for(int iz=0;iz<=res;++iz) {
    for(int iy=0;iy<=res;++iy) {
    for(int ix=0;ix<=res;++ix) {
        grid[gidx(ix,iy,iz)]=sdf(x0+ix*dx,y0+iy*dy,z0+iz*dz);
    }}}

    std::vector<V3> tris; tris.reserve(res*res*res*3);
    for(int iz=0;iz<res;++iz) {
    for(int iy=0;iy<res;++iy) {
    for(int ix=0;ix<res;++ix) {
        double gx=x0+ix*dx,gy=y0+iy*dy,gz=z0+iz*dz;
        double vv[8]; V3 vp[8];
        for(int i=0;i<8;++i){
            vv[i]=grid[gidx(ix+(int)MC_OFF[i][0],iy+(int)MC_OFF[i][1],iz+(int)MC_OFF[i][2])];
            vp[i]={gx+MC_OFF[i][0]*dx,gy+MC_OFF[i][1]*dy,gz+MC_OFF[i][2]*dz};
        }
        int ci=0;
        for(int i=0;i<8;++i) if(vv[i]<0) ci|=(1<<i);
        if(!ci||ci==255||!MC_EDGE[ci]) continue;
        V3 ep[12];
        for(int e=0;e<12;++e) if(MC_EDGE[ci]&(1<<e)) ep[e]=edge_interp(vp,vv,e);
        for(int t=0;MC_TRIS[ci][t]!=-1;t+=3){
            tris.push_back(ep[MC_TRIS[ci][t  ]]);
            tris.push_back(ep[MC_TRIS[ci][t+1]]);
            tris.push_back(ep[MC_TRIS[ci][t+2]]);
        }
    }}}
    return tris;
}

static std::string marching_cubes(
    std::function<double(double,double,double)> sdf,
    double x0,double y0,double z0,double x1,double y1,double z1,int res)
{
    return make_stl(marching_cubes_verts(sdf,x0,y0,z0,x1,y1,z1,res));
}

// ─── SDF primitives ───────────────────────────────────────────────────────────

static double sdf_box(double x,double y,double z,
                      double cx,double cy,double cz,
                      double hx,double hy,double hz){
    double dx=std::abs(x-cx)-hx, dy=std::abs(y-cy)-hy, dz=std::abs(z-cz)-hz;
    return std::min(std::max({dx,dy,dz}),0.0)
         + std::sqrt(std::max(dx,0.0)*std::max(dx,0.0)
                   + std::max(dy,0.0)*std::max(dy,0.0)
                   + std::max(dz,0.0)*std::max(dz,0.0));
}

static double sdf_ellipsoid(double x,double y,double z,
                             double cx,double cy,double cz,
                             double rx,double ry,double rz){
    double ux=(x-cx)/rx, uy=(y-cy)/ry, uz=(z-cz)/rz;
    double k0=std::sqrt(ux*ux+uy*uy+uz*uz);
    double k1=std::sqrt(ux*ux/(rx*rx)+uy*uy/(ry*ry)+uz*uz/(rz*rz));
    if(k1<1e-12) return k0-1.0;
    return k0*(k0-1.0)/k1;
}

// ─── Vertex snapping ──────────────────────────────────────────────────────────
// Project any vertex within `tol` of a box face plane exactly onto that plane.
// Vertices must also lie within the face's own bounds (± tol slack) to avoid
// false snapping on the sphere-crater surface.
static void snap_to_box_planes(std::vector<V3>& verts,
                                double cx, double cy, double cz,
                                double hx, double hy, double hz,
                                double tol)
{
    const double xlo=cx-hx, xhi=cx+hx;
    const double ylo=cy-hy, yhi=cy+hy;
    const double zlo=cz-hz, zhi=cz+hz;
    const double sl = tol; // slack for in-face bounds check

    for (auto& v : verts) {
        // ±X faces
        if (std::abs(v[0]-xlo)<tol && v[1]>=ylo-sl && v[1]<=yhi+sl && v[2]>=zlo-sl && v[2]<=zhi+sl) v[0]=xlo;
        if (std::abs(v[0]-xhi)<tol && v[1]>=ylo-sl && v[1]<=yhi+sl && v[2]>=zlo-sl && v[2]<=zhi+sl) v[0]=xhi;
        // ±Y faces
        if (std::abs(v[1]-ylo)<tol && v[0]>=xlo-sl && v[0]<=xhi+sl && v[2]>=zlo-sl && v[2]<=zhi+sl) v[1]=ylo;
        if (std::abs(v[1]-yhi)<tol && v[0]>=xlo-sl && v[0]<=xhi+sl && v[2]>=zlo-sl && v[2]<=zhi+sl) v[1]=yhi;
        // ±Z faces
        if (std::abs(v[2]-zlo)<tol && v[0]>=xlo-sl && v[0]<=xhi+sl && v[1]>=ylo-sl && v[1]<=yhi+sl) v[2]=zlo;
        if (std::abs(v[2]-zhi)<tol && v[0]>=xlo-sl && v[0]<=xhi+sl && v[1]>=ylo-sl && v[1]<=yhi+sl) v[2]=zhi;
    }
}

// ─── Ellipsoid surface projection ─────────────────────────────────────────────
// After snap+smooth, marching-cubes verts on the sphere crater are still slightly
// off the true ellipsoid surface. This projects each such vertex exactly onto the
// ellipsoid by scaling its direction from the ellipsoid centre.
// Vertices near any box face plane are left untouched (they are already snapped
// flat and define the sharp box–sphere boundary edge).
static void project_to_ellipsoid(std::vector<V3>& tris,
    double bcx,double bcy,double bcz,
    double bhx,double bhy,double bhz,
    double scx,double scy,double scz,
    double srx,double sry,double srz,
    double face_margin,   // protect verts this close to any box face plane
    double sph_tol)       // project verts whose |ellipsoid_sdf| < this
{
    const double xlo=bcx-bhx, xhi=bcx+bhx;
    const double ylo=bcy-bhy, yhi=bcy+bhy;
    const double zlo=bcz-bhz, zhi=bcz+bhz;

    for (auto& v : tris) {
        // Leave vertices near any box face plane alone
        if (std::abs(v[0]-xlo)<face_margin || std::abs(v[0]-xhi)<face_margin) continue;
        if (std::abs(v[1]-ylo)<face_margin || std::abs(v[1]-yhi)<face_margin) continue;
        if (std::abs(v[2]-zlo)<face_margin || std::abs(v[2]-zhi)<face_margin) continue;

        // Only project if close to the ellipsoid surface
        double es = sdf_ellipsoid(v[0],v[1],v[2], scx,scy,scz, srx,sry,srz);
        if (std::abs(es) > sph_tol) continue;

        // Scale direction from ellipsoid centre so the point lands on the surface
        double ux=(v[0]-scx)/srx, uy=(v[1]-scy)/sry, uz=(v[2]-scz)/srz;
        double len=std::sqrt(ux*ux+uy*uy+uz*uz);
        if (len<1e-12) continue;
        v[0]=scx+(ux/len)*srx;
        v[1]=scy+(uy/len)*sry;
        v[2]=scz+(uz/len)*srz;
    }
}

// Helper: drop degenerate/NaN triangles in-place
static void filter_degenerate(std::vector<V3>& tris) {
    std::vector<V3> clean; clean.reserve(tris.size());
    auto fin3=[](const V3& v){
        return std::isfinite(v[0])&&std::isfinite(v[1])&&std::isfinite(v[2]);
    };
    for(int i=0;i<(int)tris.size();i+=3){
        const V3& A=tris[i],&B=tris[i+1],&C=tris[i+2];
        if(!fin3(A)||!fin3(B)||!fin3(C)) continue;
        V3 ab={B[0]-A[0],B[1]-A[1],B[2]-A[2]};
        V3 ac={C[0]-A[0],C[1]-A[1],C[2]-A[2]};
        V3 cr={ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]};
        if(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]>1e-20)
            { clean.push_back(A); clean.push_back(B); clean.push_back(C); }
    }
    tris=std::move(clean);
}

// ─── Boolean subtract ─────────────────────────────────────────────────────────

std::string boolean_subtract(double bcx,double bcy,double bcz,
                              double bhx,double bhy,double bhz,
                              double scx,double scy,double scz,
                              double srx,double sry,double srz,
                              int resolution)
{
    double xlo=std::min(bcx-bhx,scx-srx)-0.05;
    double ylo=std::min(bcy-bhy,scy-sry)-0.05;
    double zlo=std::min(bcz-bhz,scz-srz)-0.05;
    double xhi=std::max(bcx+bhx,scx+srx)+0.05;
    double yhi=std::max(bcy+bhy,scy+sry)+0.05;
    double zhi=std::max(bcz+bhz,scz+srz)+0.05;

    auto sdf=[=](double x,double y,double z){
        double a=sdf_box      (x,y,z, bcx,bcy,bcz, bhx,bhy,bhz);
        double b=sdf_ellipsoid(x,y,z, scx,scy,scz, srx,sry,srz);
        return std::max(a,-b);
    };

    std::vector<V3> tris=marching_cubes_verts(sdf,xlo,ylo,zlo,xhi,yhi,zhi,resolution);

    // Grid step — MC verts can be at most 0.5 steps off a true planar surface
    double step = std::max({(xhi-xlo)/resolution,
                            (yhi-ylo)/resolution,
                            (zhi-zlo)/resolution});

    // Pass 1: snap near-plane vertices to exact box face planes
    snap_to_box_planes(tris, bcx,bcy,bcz, bhx,bhy,bhz, step * 0.52);

    // More iterations + higher alpha to push MC staircase vertices closer to the
    // true surface before the analytic projection step below.
    smooth_mesh(tris, 20, 0.70, 30.0);

    // Pass 2: tighter re-snap corrects any smoothing drift on flat faces
    snap_to_box_planes(tris, bcx,bcy,bcz, bhx,bhy,bhz, step * 0.12);

    // KEY FIX: project every sphere-crater vertex exactly onto the analytic ellipsoid.
    // After this step the sphere surface is geometrically perfect — the cross-product
    // normals written by make_stl are then analytically correct with no staircase error.
    // Vertices within face_margin of any box plane are protected so the sharp
    // box–sphere boundary edge is preserved.
    project_to_ellipsoid(tris, bcx,bcy,bcz, bhx,bhy,bhz,
                               scx,scy,scz, srx,sry,srz,
                               step * 0.35,   // face_margin
                               step * 2.0);   // sph_tol

    // Remove degenerate triangles (projection can make near-coincident verts identical)
    filter_degenerate(tris);

    return make_stl(tris);
}

// ─── Cylinder ─────────────────────────────────────────────────────────────────

std::string cylinder_stl(double cx, double cy, double cz,
                          double radius, double height, int segs)
{
    const double pi = std::acos(-1.0);
    std::vector<V3> t;
    t.reserve(segs * 4 * 3);

    double ytop = cy + height * 0.5;
    double ybot = cy - height * 0.5;
    V3 top_center  = {cx, ytop, cz};
    V3 bot_center  = {cx, ybot, cz};

    for (int i = 0; i < segs; ++i) {
        double a0 = 2.0 * pi * i       / segs;
        double a1 = 2.0 * pi * (i + 1) / segs;
        double c0 = std::cos(a0), s0 = std::sin(a0);
        double c1 = std::cos(a1), s1 = std::sin(a1);

        V3 t0 = {cx + radius * c0, ytop, cz + radius * s0};
        V3 t1 = {cx + radius * c1, ytop, cz + radius * s1};
        V3 b0 = {cx + radius * c0, ybot, cz + radius * s0};
        V3 b1 = {cx + radius * c1, ybot, cz + radius * s1};

        // Top cap (normal +Y) — wind CCW when viewed from above
        t.push_back(top_center); t.push_back(t1); t.push_back(t0);
        // Bottom cap (normal -Y) — wind CCW when viewed from below
        t.push_back(bot_center); t.push_back(b0); t.push_back(b1);
        // Side quad: two triangles, outward radial normals
        t.push_back(b0); t.push_back(t0); t.push_back(t1);
        t.push_back(b0); t.push_back(t1); t.push_back(b1);
    }
    return make_stl(t);
}

// ─── Cone ─────────────────────────────────────────────────────────────────────

std::string cone_stl(double cx, double cy, double cz,
                     double radius, double height, int segs)
{
    const double pi = std::acos(-1.0);
    std::vector<V3> t;
    t.reserve(segs * 2 * 3);

    double ytop = cy + height * 0.5;  // apex
    double ybot = cy - height * 0.5;  // base center
    V3 apex       = {cx, ytop, cz};
    V3 bot_center = {cx, ybot, cz};

    for (int i = 0; i < segs; ++i) {
        double a0 = 2.0 * pi * i       / segs;
        double a1 = 2.0 * pi * (i + 1) / segs;
        double c0 = std::cos(a0), s0 = std::sin(a0);
        double c1 = std::cos(a1), s1 = std::sin(a1);

        V3 b0 = {cx + radius * c0, ybot, cz + radius * s0};
        V3 b1 = {cx + radius * c1, ybot, cz + radius * s1};

        // Side triangle: base edge to apex (outward-slanted normal via cross product)
        t.push_back(b0); t.push_back(b1); t.push_back(apex);
        // Bottom cap fan (normal -Y): wind so normal points downward
        t.push_back(bot_center); t.push_back(b1); t.push_back(b0);
    }
    return make_stl(t);
}

// ─── Torus ────────────────────────────────────────────────────────────────────

std::string torus_stl(double cx, double cy, double cz,
                       double major_r, double minor_r,
                       int major_segs, int minor_segs)
{
    const double pi = std::acos(-1.0);
    std::vector<V3> t;
    t.reserve(major_segs * minor_segs * 2 * 3);

    // Parameterization: u around the ring (major), v around the tube (minor)
    // Ring in XZ plane, Y is up
    auto point = [&](int ui, int vi) -> V3 {
        double u = 2.0 * pi * ui / major_segs;
        double v = 2.0 * pi * vi / minor_segs;
        double cu = std::cos(u), su = std::sin(u);
        double cv = std::cos(v), sv = std::sin(v);
        // Tube center is at (major_r*cu, 0, major_r*su) in XZ plane (Y=0 ring)
        double x = cx + (major_r + minor_r * cv) * cu;
        double y = cy + minor_r * sv;
        double z = cz + (major_r + minor_r * cv) * su;
        return {x, y, z};
    };

    for (int ui = 0; ui < major_segs; ++ui) {
        for (int vi = 0; vi < minor_segs; ++vi) {
            V3 p00 = point(ui,     vi    );
            V3 p10 = point(ui + 1, vi    );
            V3 p01 = point(ui,     vi + 1);
            V3 p11 = point(ui + 1, vi + 1);
            // Two triangles per quad, outward normals computed by write_tri
            t.push_back(p00); t.push_back(p10); t.push_back(p11);
            t.push_back(p00); t.push_back(p11); t.push_back(p01);
        }
    }
    return make_stl(t);
}

// ─── Plane ────────────────────────────────────────────────────────────────────

std::string plane_stl(double cx, double cy, double cz,
                       double width, double depth, int subdiv)
{
    std::vector<V3> t;
    t.reserve(subdiv * subdiv * 2 * 3);

    double half_w = width * 0.5;
    double half_d = depth * 0.5;

    for (int i = 0; i < subdiv; ++i) {
        for (int j = 0; j < subdiv; ++j) {
            double x0 = cx - half_w + width  * i       / subdiv;
            double x1 = cx - half_w + width  * (i + 1) / subdiv;
            double z0 = cz - half_d + depth  * j       / subdiv;
            double z1 = cz - half_d + depth  * (j + 1) / subdiv;

            V3 v00 = {x0, cy, z0};
            V3 v10 = {x1, cy, z0};
            V3 v01 = {x0, cy, z1};
            V3 v11 = {x1, cy, z1};

            // Wind CCW when viewed from above (+Y normal)
            t.push_back(v00); t.push_back(v10); t.push_back(v11);
            t.push_back(v00); t.push_back(v11); t.push_back(v01);
        }
    }
    return make_stl(t);
}

// ─── Icosphere ───────────────────────────────────────────────────────────────
// Subdivided icosahedron, all triangles roughly equal size.

std::string icosphere_stl(double cx, double cy, double cz,
                           double radius, int subdivisions)
{
    const double X = 0.525731112119133606;
    const double Z = 0.850650808352039932;
    // 12 icosahedron vertices on unit sphere
    std::vector<V3> verts = {
        {-X,0,Z},{X,0,Z},{-X,0,-Z},{X,0,-Z},
        {0,Z,X},{0,Z,-X},{0,-Z,X},{0,-Z,-X},
        {Z,X,0},{-Z,X,0},{Z,-X,0},{-Z,-X,0}
    };
    // 20 icosahedron faces
    std::vector<std::array<int,3>> faces = {
        {0,4,1},{0,9,4},{9,5,4},{4,5,8},{4,8,1},
        {8,10,1},{8,3,10},{5,3,8},{5,2,3},{2,7,3},
        {7,10,3},{7,6,10},{7,11,6},{11,0,6},{0,1,6},
        {6,1,10},{9,0,11},{9,11,2},{9,2,5},{7,2,11}
    };

    // Midpoint cache to avoid duplicates
    std::map<std::pair<int,int>, int> midCache;
    auto midpoint = [&](int a, int b) -> int {
        auto key = std::make_pair(std::min(a,b), std::max(a,b));
        auto it = midCache.find(key);
        if (it != midCache.end()) return it->second;
        V3 va = verts[a], vb = verts[b];
        V3 m = {(va[0]+vb[0])*0.5, (va[1]+vb[1])*0.5, (va[2]+vb[2])*0.5};
        double len = std::sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
        if (len > 1e-12) { m[0]/=len; m[1]/=len; m[2]/=len; }
        int idx = (int)verts.size();
        verts.push_back(m);
        midCache[key] = idx;
        return idx;
    };

    for (int s = 0; s < subdivisions; ++s) {
        std::vector<std::array<int,3>> next;
        next.reserve(faces.size() * 4);
        for (auto& f : faces) {
            int a = midpoint(f[0], f[1]);
            int b = midpoint(f[1], f[2]);
            int c = midpoint(f[2], f[0]);
            next.push_back({f[0], a, c});
            next.push_back({f[1], b, a});
            next.push_back({f[2], c, b});
            next.push_back({a,    b, c});
        }
        faces = std::move(next);
    }

    // Scale and translate, then emit triangles
    std::vector<V3> t;
    t.reserve(faces.size() * 3);
    for (auto& f : faces) {
        V3 a = {cx + verts[f[0]][0]*radius, cy + verts[f[0]][1]*radius, cz + verts[f[0]][2]*radius};
        V3 b = {cx + verts[f[1]][0]*radius, cy + verts[f[1]][1]*radius, cz + verts[f[1]][2]*radius};
        V3 c = {cx + verts[f[2]][0]*radius, cy + verts[f[2]][1]*radius, cz + verts[f[2]][2]*radius};
        t.push_back(a); t.push_back(b); t.push_back(c);
    }
    return make_stl(t);
}

// ---------------------------------------------------------------------------
// Public wrappers — expose MC internals to other translation units
// ---------------------------------------------------------------------------

std::vector<MCV3> marching_cubes_run(
    std::function<double(double,double,double)> sdf,
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    int res)
{
    return marching_cubes_verts(std::move(sdf), x0, y0, z0, x1, y1, z1, res);
}

SurfaceManifold tris_to_manifold(const std::vector<MCV3>& tris)
{
    SurfaceManifold m;
    const uint32_t n = static_cast<uint32_t>(tris.size() / 3);
    m.vertices.reserve(tris.size());
    m.faces.reserve(n);

    for (uint32_t i = 0; i < n; ++i) {
        const MCV3& A = tris[i*3+0];
        const MCV3& B = tris[i*3+1];
        const MCV3& C = tris[i*3+2];

        uint32_t base = static_cast<uint32_t>(m.vertices.size());
        m.vertices.push_back({{A[0], A[1], A[2]}});
        m.vertices.push_back({{B[0], B[1], B[2]}});
        m.vertices.push_back({{C[0], C[1], C[2]}});

        // Face normal
        double ux=B[0]-A[0], uy=B[1]-A[1], uz=B[2]-A[2];
        double vx=C[0]-A[0], vy=C[1]-A[1], vz=C[2]-A[2];
        double nx=uy*vz-uz*vy, ny=uz*vx-ux*vz, nz=ux*vy-uy*vx;
        double len=std::sqrt(nx*nx+ny*ny+nz*nz);
        if(len > 1e-12){ nx/=len; ny/=len; nz/=len; }
        else { nx=0; ny=0; nz=1; }

        Triangle tri;
        tri.vi     = {base, base+1, base+2};
        tri.normal = {nx, ny, nz};
        m.faces.push_back(tri);
    }
    return m;
}

} // namespace cfd
