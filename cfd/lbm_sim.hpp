#pragma once
#include <vector>
#include <cstdint>

struct SimObstacle {
    float cx, cy, r; // normalized 0..1 coords
};

struct SimParams {
    int   nx        = 200;
    int   ny        = 80;
    float reynolds  = 150.f;
    float inlet_vel = 0.1f;
    int   steps     = 4000;
    std::vector<SimObstacle> obstacles;
};

struct SimResult {
    int nx, ny, steps_done;
    std::vector<float>   vx, vy, speed, rho;
    std::vector<uint8_t> solid;
};

SimResult run_lbm(const SimParams& p);
