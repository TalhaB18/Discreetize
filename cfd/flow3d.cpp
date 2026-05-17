// 3-D potential-flow streamline tracer
// Uniform inflow + sphere-doublet perturbation per obstacle, RK4 integration
#include "flow3d.hpp"
#include <cmath>
#include <algorithm>
#include <cstring>

struct V3 { float x,y,z; };
static V3 add(V3 a,V3 b){return{a.x+b.x,a.y+b.y,a.z+b.z};}
static V3 scl(V3 a,float s){return{a.x*s,a.y*s,a.z*s};}
static float dot(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static float len(V3 a){return std::sqrt(dot(a,a));}

struct Sphere { float cx,cy,cz,R; };

// Potential flow velocity: uniform + doublet perturbation from each sphere
// Reference: Lamb, Hydrodynamics §100
static V3 vel_at(V3 p, float U, float fd[3], const std::vector<Sphere>& obs) {
    V3 v = {U*fd[0], U*fd[1], U*fd[2]};
    for (const auto& s : obs) {
        V3 d = {p.x-s.cx, p.y-s.cy, p.z-s.cz};
        float r2 = dot(d,d);
        float r  = std::sqrt(r2);
        if (r < s.R * 0.05f) return {0,0,0}; // inside / stagnation
        float R3  = s.R*s.R*s.R;
        float r3  = r*r2, r5 = r3*r2;
        // dot(flow_dir, d)
        float dfd = fd[0]*d.x + fd[1]*d.y + fd[2]*d.z;
        // Doublet: phi = (U*R3/2) * dfd / r3
        // grad phi_pert:
        float f = U * R3 * 0.5f;
        v.x += f * (fd[0]/r3 - 3.f*dfd*d.x/r5);
        v.y += f * (fd[1]/r3 - 3.f*dfd*d.y/r5);
        v.z += f * (fd[2]/r3 - 3.f*dfd*d.z/r5);
    }
    return v;
}

static V3 rk4_step(V3 p, float h, float U, float fd[3], const std::vector<Sphere>& obs) {
    V3 k1 = vel_at(p,                  U, fd, obs);
    V3 k2 = vel_at(add(p,scl(k1,h*.5f)),U, fd, obs);
    V3 k3 = vel_at(add(p,scl(k2,h*.5f)),U, fd, obs);
    V3 k4 = vel_at(add(p,scl(k3,h)),   U, fd, obs);
    V3 dv = add(add(k1,k4), scl(add(k2,k3),2.f));
    return add(p, scl(dv, h/6.f));
}

FlowResult compute_flow(const FlowParams& p) {
    // Build equivalent spheres from boxes
    std::vector<Sphere> obs;
    obs.reserve(p.objects.size());
    for (const auto& b : p.objects) {
        float vol = b.sx * b.sy * b.sz;
        float R   = 0.5f * std::cbrt(vol);
        obs.push_back({b.cx, b.cy, b.cz, std::max(R, 0.05f)});
    }

    // Domain bounds — encompass all objects + generous padding
    float mn[3]={1e9,1e9,1e9}, mx[3]={-1e9,-1e9,-1e9};
    if (p.objects.empty()) {
        float d[3]={-5,-3,-3}, u[3]={8,3,3};
        memcpy(mn,d,12); memcpy(mx,u,12);
    } else {
        for (const auto& b : p.objects) {
            mn[0]=std::min(mn[0],b.cx-b.sx*.5f);
            mn[1]=std::min(mn[1],b.cy-b.sy*.5f);
            mn[2]=std::min(mn[2],b.cz-b.sz*.5f);
            mx[0]=std::max(mx[0],b.cx+b.sx*.5f);
            mx[1]=std::max(mx[1],b.cy+b.sy*.5f);
            mx[2]=std::max(mx[2],b.cz+b.sz*.5f);
        }
        float sz[3]={mx[0]-mn[0],mx[1]-mn[1],mx[2]-mn[2]};
        float maxS = std::max({sz[0],sz[1],sz[2],0.5f});
        // Long upstream/downstream so streams visibly start far from objects
        mn[0]-=maxS*3.0f; mx[0]+=maxS*3.0f;
        mn[1]-=maxS*1.4f; mx[1]+=maxS*1.4f;
        mn[2]-=maxS*1.4f; mx[2]+=maxS*1.4f;
    }

    float ds[3]={mx[0]-mn[0],mx[1]-mn[1],mx[2]-mn[2]};
    float domLen = std::max({ds[0],ds[1],ds[2]});
    // Step size: cross the domain in max_steps
    float h = domLen / (float)p.max_steps;

    // Seed grid perpendicular to flow direction on the inlet face
    // For fd=(1,0,0): seed on x=mn[0]+5%, grid over y and z
    int nsq = (int)std::ceil(std::sqrt((double)p.seed_count));

    FlowResult res;
    res.vmax = 0;
    memcpy(res.domain,   mn, 12);
    memcpy(res.domain+3, mx, 12);

    float fd[3];
    memcpy(fd, p.flow_dir, 12);
    // Normalise flow dir
    float flen = std::sqrt(fd[0]*fd[0]+fd[1]*fd[1]+fd[2]*fd[2]);
    if (flen < 1e-6f) { fd[0]=1; fd[1]=fd[2]=0; flen=1; }
    fd[0]/=flen; fd[1]/=flen; fd[2]/=flen;

    // Seed axes (two axes perpendicular to flow dir)
    // Simple cross with world-up or world-right
    V3 fv = {fd[0],fd[1],fd[2]};
    V3 up = (std::abs(fd[1]) < 0.9f) ? V3{0,1,0} : V3{1,0,0};
    // right = fv × up
    V3 right = {fv.y*up.z - fv.z*up.y, fv.z*up.x - fv.x*up.z, fv.x*up.y - fv.y*up.x};
    float rlen = len(right); if (rlen<1e-6f) rlen=1;
    right = scl(right,1.f/rlen);
    // up2 = right × fv
    V3 up2 = {right.y*fv.z - right.z*fv.y, right.z*fv.x - right.x*fv.z, right.x*fv.y - right.y*fv.x};

    // Center of inlet face: mn + fd * 2% of domain
    float inletDist = 0.02f * domLen;
    float cx = (mn[0]+mx[0])*.5f - fd[0]*ds[0]*.45f + fd[0]*inletDist;
    float cy = (mn[1]+mx[1])*.5f - fd[1]*ds[1]*.45f + fd[1]*inletDist;
    float cz = (mn[2]+mx[2])*.5f - fd[2]*ds[2]*.45f + fd[2]*inletDist;

    float span1 = std::max(ds[1], ds[2]) * 0.85f;
    float span2 = span1;

    for (int i = 0; i < nsq && (int)res.lines.size() < p.seed_count; ++i) {
        for (int j = 0; j < nsq && (int)res.lines.size() < p.seed_count; ++j) {
            float u1 = (nsq>1) ? (float)i/(nsq-1)-0.5f : 0.f;
            float u2 = (nsq>1) ? (float)j/(nsq-1)-0.5f : 0.f;
            V3 seed = {
                cx + right.x*u1*span1 + up2.x*u2*span2,
                cy + right.y*u1*span1 + up2.y*u2*span2,
                cz + right.z*u1*span1 + up2.z*u2*span2
            };

            StreamLine line;
            V3 cur = seed;
            for (int s = 0; s < p.max_steps; ++s) {
                if (cur.x<mn[0]||cur.x>mx[0]||
                    cur.y<mn[1]||cur.y>mx[1]||
                    cur.z<mn[2]||cur.z>mx[2]) break;
                V3 v = vel_at(cur, p.inlet_speed, fd, obs);
                float sp = len(v);
                if (sp < 1e-6f) break;
                line.pts.push_back(cur.x);
                line.pts.push_back(cur.y);
                line.pts.push_back(cur.z);
                line.speed.push_back(sp);
                if (sp > res.vmax) res.vmax = sp;
                // Adaptive step: advance fixed spatial distance
                float ht = h / sp;
                cur = rk4_step(cur, ht, p.inlet_speed, fd, obs);
            }
            if (line.pts.size() >= 9)
                res.lines.push_back(std::move(line));
        }
    }

    if (res.vmax < 1e-9f) res.vmax = p.inlet_speed;
    return res;
}
