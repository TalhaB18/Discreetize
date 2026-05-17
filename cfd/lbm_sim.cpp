// D2Q9 Lattice-Boltzmann flow solver (BGK, bounce-back walls, Zou-He inlet)
#include "lbm_sim.hpp"
#include <cmath>
#include <algorithm>

static constexpr int QQ = 9;
static constexpr int EX[9] = { 0, 1, 0,-1, 0, 1,-1,-1, 1 };
static constexpr int EY[9] = { 0, 0, 1, 0,-1, 1, 1,-1,-1 };
static constexpr int OP[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };
static constexpr float WW[9] = {
    4.f/9.f,
    1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f,
    1.f/36.f,1.f/36.f,1.f/36.f,1.f/36.f
};

static inline float feq(int q, float rho, float ux, float uy) {
    float eu = EX[q]*ux + EY[q]*uy;
    return WW[q] * rho * (1.f + 3.f*eu + 4.5f*eu*eu - 1.5f*(ux*ux + uy*uy));
}

SimResult run_lbm(const SimParams& p) {
    const int NX = p.nx, NY = p.ny, N = NX * NY;

    // ── Solid mask ───────────────────────────────────────────────
    std::vector<uint8_t> solid(N, 0);
    for (int x = 0; x < NX; ++x) {
        solid[0*NX + x]       = 1; // bottom wall
        solid[(NY-1)*NX + x]  = 1; // top wall
    }
    for (auto& o : p.obstacles) {
        int ocx = (int)(o.cx * NX);
        int ocy = (int)(o.cy * NY);
        int orr = std::max(2, (int)(o.r * std::min(NX, NY)));
        for (int y = 0; y < NY; ++y)
            for (int x = 0; x < NX; ++x)
                if ((x-ocx)*(x-ocx)+(y-ocy)*(y-ocy) <= orr*orr)
                    solid[y*NX+x] = 1;
    }

    // ── Relaxation ───────────────────────────────────────────────
    float nu    = p.inlet_vel * (float)NY / p.reynolds;
    float tau   = std::max(0.52f, 3.f*nu + 0.5f);
    float omega = 1.f / tau;

    // ── Distribution functions (flat: f[cell*QQ + q]) ────────────
    std::vector<float> f(N*QQ), g(N*QQ);
    for (int i = 0; i < N; ++i)
        if (!solid[i])
            for (int q = 0; q < QQ; ++q)
                f[i*QQ+q] = WW[q]; // rest equilibrium (rho=1, u=0)

    // ── Main loop ────────────────────────────────────────────────
    for (int step = 0; step < p.steps; ++step) {

        // 1. Streaming (pull scheme)
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                int i = y*NX + x;
                for (int q = 0; q < QQ; ++q) {
                    int sx = x - EX[q], sy = y - EY[q];
                    if (sx < 0 || sx >= NX || sy < 0 || sy >= NY) {
                        g[i*QQ+q] = WW[q]; // open boundary default
                        continue;
                    }
                    int si = sy*NX + sx;
                    // Bounce-back: if source is solid, reflect
                    g[i*QQ+q] = solid[si] ? f[i*QQ + OP[q]] : f[si*QQ + q];
                }
            }
        }
        f.swap(g);

        // 2. Inlet BC (x=0): Zou-He fixed velocity
        float ux0 = p.inlet_vel, uy0 = 0.f;
        for (int y = 1; y < NY-1; ++y) {
            if (solid[y*NX]) continue;
            float* fi = &f[y*NX*QQ];
            float rho = (fi[0]+fi[2]+fi[4] + 2.f*(fi[3]+fi[6]+fi[7])) / (1.f - ux0);
            for (int q = 0; q < QQ; ++q) fi[q] = feq(q, rho, ux0, uy0);
        }

        // 3. Outlet BC (x=NX-1): zero-gradient
        for (int y = 1; y < NY-1; ++y) {
            float* dst = &f[(y*NX + NX-1)*QQ];
            float* src = &f[(y*NX + NX-2)*QQ];
            for (int q = 0; q < QQ; ++q) dst[q] = src[q];
        }

        // 4. BGK collision
        for (int i = 0; i < N; ++i) {
            if (solid[i]) continue;
            float rho = 0, ux = 0, uy = 0;
            for (int q = 0; q < QQ; ++q) {
                rho += f[i*QQ+q];
                ux  += EX[q]*f[i*QQ+q];
                uy  += EY[q]*f[i*QQ+q];
            }
            ux /= rho; uy /= rho;
            for (int q = 0; q < QQ; ++q)
                f[i*QQ+q] += omega * (feq(q, rho, ux, uy) - f[i*QQ+q]);
        }
    }

    // ── Extract result ───────────────────────────────────────────
    SimResult res;
    res.nx = NX; res.ny = NY; res.steps_done = p.steps;
    res.vx.resize(N); res.vy.resize(N);
    res.speed.resize(N); res.rho.resize(N);
    res.solid.assign(solid.begin(), solid.end());

    for (int i = 0; i < N; ++i) {
        if (solid[i]) { res.vx[i]=res.vy[i]=res.speed[i]=res.rho[i]=0; continue; }
        float rho=0, ux=0, uy=0;
        for (int q = 0; q < QQ; ++q) {
            rho += f[i*QQ+q];
            ux  += EX[q]*f[i*QQ+q];
            uy  += EY[q]*f[i*QQ+q];
        }
        ux /= rho; uy /= rho;
        res.rho[i]   = rho;
        res.vx[i]    = ux;
        res.vy[i]    = uy;
        res.speed[i] = std::sqrt(ux*ux + uy*uy);
    }
    return res;
}
