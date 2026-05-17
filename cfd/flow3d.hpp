#pragma once
#include <vector>

struct FlowBox {
    float cx, cy, cz;  // center
    float sx, sy, sz;  // full extents
};

struct StreamLine {
    std::vector<float> pts;   // x,y,z,x,y,z,... flat
    std::vector<float> speed; // speed per point
};

struct FlowResult {
    std::vector<StreamLine> lines;
    float vmax;
    float domain[6]; // xmin,ymin,zmin,xmax,ymax,zmax
};

struct FlowParams {
    std::vector<FlowBox> objects;
    float flow_dir[3]  = {1, 0, 0};
    float inlet_speed  = 1.0f;
    int   seed_count   = 200;
    int   max_steps    = 400;
};

FlowResult compute_flow(const FlowParams& p);
