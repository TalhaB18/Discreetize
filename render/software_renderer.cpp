// =============================================================================
// Software rasterizer with Phong shading → PNG output
// No external image library dependencies.
// =============================================================================

#include "software_renderer.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace render {

// =============================================================================
// Math helpers
// =============================================================================

struct Vec3 {
    double x, y, z;
    Vec3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& o) const { return {x+o.x, y+o.y, z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x, y-o.y, z-o.z}; }
    Vec3 operator*(double s)      const { return {x*s,   y*s,   z*s};   }
    Vec3 operator/(double s)      const { return {x/s,   y/s,   z/s};   }
    double dot(const Vec3& o)     const { return x*o.x + y*o.y + z*o.z; }
    Vec3 cross(const Vec3& o)     const {
        return {y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x};
    }
    double len2() const { return dot(*this); }
    double len()  const { return std::sqrt(len2()); }
    Vec3 norm()   const { double l = len(); return (l > 1e-12) ? (*this / l) : Vec3(0,0,0); }
};

// 4×4 matrix, row-major
struct Mat4 {
    double m[4][4] = {};
    static Mat4 identity() {
        Mat4 r; r.m[0][0]=r.m[1][1]=r.m[2][2]=r.m[3][3]=1; return r;
    }
    Mat4 operator*(const Mat4& o) const {
        Mat4 r;
        for (int i=0;i<4;++i)
            for (int j=0;j<4;++j)
                for (int k=0;k<4;++k)
                    r.m[i][j] += m[i][k]*o.m[k][j];
        return r;
    }
    // Multiply homogeneous point [x,y,z,1] → returns Vec3 after w-divide
    Vec3 mulPoint(const Vec3& v) const {
        double x = m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z + m[0][3];
        double y = m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z + m[1][3];
        double z = m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z + m[2][3];
        double w = m[3][0]*v.x + m[3][1]*v.y + m[3][2]*v.z + m[3][3];
        if (std::fabs(w) > 1e-12) { x/=w; y/=w; z/=w; }
        return {x, y, z};
    }
    // Multiply direction (w=0)
    Vec3 mulDir(const Vec3& v) const {
        double x = m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z;
        double y = m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z;
        double z = m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z;
        return {x, y, z};
    }
};

// Build look-at view matrix (camera-space transform)
static Mat4 look_at(Vec3 eye, Vec3 target, Vec3 up) {
    Vec3 f = (target - eye).norm();   // forward
    Vec3 r = f.cross(up).norm();      // right
    Vec3 u = r.cross(f);              // up (recomputed)

    Mat4 M;
    M.m[0][0]=r.x;  M.m[0][1]=r.y;  M.m[0][2]=r.z;  M.m[0][3]=-r.dot(eye);
    M.m[1][0]=u.x;  M.m[1][1]=u.y;  M.m[1][2]=u.z;  M.m[1][3]=-u.dot(eye);
    M.m[2][0]=-f.x; M.m[2][1]=-f.y; M.m[2][2]=-f.z; M.m[2][3]=f.dot(eye);
    M.m[3][3]=1;
    return M;
}

// Build perspective projection matrix (OpenGL convention, NDC -1..1)
static Mat4 perspective(double fov_rad, double aspect, double near_z, double far_z) {
    double t = std::tan(fov_rad * 0.5);
    Mat4 M;
    M.m[0][0] = 1.0 / (aspect * t);
    M.m[1][1] = 1.0 / t;
    M.m[2][2] = -(far_z + near_z) / (far_z - near_z);
    M.m[2][3] = -(2.0 * far_z * near_z) / (far_z - near_z);
    M.m[3][2] = -1.0;
    return M;
}

// =============================================================================
// PNG encoder — no external dependencies
// =============================================================================

// CRC32 table
static uint32_t crc_table[256];
static bool     crc_table_ready = false;

static void build_crc_table() {
    if (crc_table_ready) return;
    for (uint32_t n = 0; n < 256; ++n) {
        uint32_t c = n;
        for (int k = 0; k < 8; ++k)
            c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
        crc_table[n] = c;
    }
    crc_table_ready = true;
}

static uint32_t crc32(const uint8_t* buf, size_t len, uint32_t crc = 0xFFFFFFFFu) {
    build_crc_table();
    for (size_t i = 0; i < len; ++i)
        crc = crc_table[(crc ^ buf[i]) & 0xFF] ^ (crc >> 8);
    return crc ^ 0xFFFFFFFFu;
}

static void push_u32be(std::vector<uint8_t>& v, uint32_t x) {
    v.push_back((x >> 24) & 0xFF);
    v.push_back((x >> 16) & 0xFF);
    v.push_back((x >>  8) & 0xFF);
    v.push_back((x      ) & 0xFF);
}

static void push_u16le(std::vector<uint8_t>& v, uint16_t x) {
    v.push_back( x       & 0xFF);
    v.push_back((x >> 8) & 0xFF);
}

static void write_chunk(std::vector<uint8_t>& out,
                        const char type[4],
                        const std::vector<uint8_t>& data) {
    push_u32be(out, static_cast<uint32_t>(data.size()));
    size_t crc_start = out.size();
    out.push_back(type[0]); out.push_back(type[1]);
    out.push_back(type[2]); out.push_back(type[3]);
    for (auto b : data) out.push_back(b);
    uint32_t c = crc32(out.data() + crc_start, 4 + data.size());
    push_u32be(out, c);
}

// Adler-32 for zlib wrapper
static uint32_t adler32(const uint8_t* buf, size_t len) {
    uint32_t s1 = 1, s2 = 0;
    for (size_t i = 0; i < len; ++i) {
        s1 = (s1 + buf[i]) % 65521u;
        s2 = (s2 + s1)     % 65521u;
    }
    return (s2 << 16) | s1;
}

// Encode raw bytes using zlib stored blocks (no compression, but valid zlib)
static std::vector<uint8_t> zlib_stored(const uint8_t* data, size_t len) {
    std::vector<uint8_t> out;
    // Zlib header: CMF=0x78, FLG such that (CMF*256+FLG) % 31 == 0
    // 0x78 * 256 + 0x01 = 30721; 30721 % 31 = 0  ✓
    out.push_back(0x78);
    out.push_back(0x01);

    // Stored blocks: each block can hold up to 65535 bytes
    const size_t BLOCK = 65535;
    size_t offset = 0;
    while (offset < len || len == 0) {
        size_t chunk = std::min(BLOCK, len - offset);
        bool   final = (offset + chunk >= len);
        out.push_back(final ? 0x01 : 0x00); // BFINAL | BTYPE=00 (no compress)
        uint16_t nlen  = static_cast<uint16_t>(chunk);
        uint16_t nlenC = static_cast<uint16_t>(~nlen);
        push_u16le(out, nlen);
        push_u16le(out, nlenC);
        for (size_t i = 0; i < chunk; ++i)
            out.push_back(data[offset + i]);
        offset += chunk;
        if (final) break;
    }

    uint32_t a = adler32(data, len);
    push_u32be(out, a); // big-endian Adler-32
    return out;
}

// Encode RGBA image to PNG bytes.
// pixels: row-major, 4 bytes per pixel (RGBA)
static std::vector<uint8_t> encode_png(const std::vector<uint8_t>& pixels,
                                       int width, int height) {
    std::vector<uint8_t> out;
    // PNG signature
    const uint8_t sig[8] = {137,80,78,71,13,10,26,10};
    for (auto b : sig) out.push_back(b);

    // IHDR
    {
        std::vector<uint8_t> ihdr(13);
        ihdr[0]=(width>>24)&0xFF; ihdr[1]=(width>>16)&0xFF;
        ihdr[2]=(width>> 8)&0xFF; ihdr[3]=(width    )&0xFF;
        ihdr[4]=(height>>24)&0xFF; ihdr[5]=(height>>16)&0xFF;
        ihdr[6]=(height>> 8)&0xFF; ihdr[7]=(height    )&0xFF;
        ihdr[8]=8;  // bit depth
        ihdr[9]=2;  // color type: RGB (not RGBA to keep it simple with filter byte)
        ihdr[10]=0; ihdr[11]=0; ihdr[12]=0;
        write_chunk(out, "IHDR", ihdr);
    }

    // IDAT — build raw scan-line data (filter byte 0 per row) then compress
    {
        // We use RGB (3 bytes per pixel, color type 2)
        std::vector<uint8_t> raw;
        raw.reserve(static_cast<size_t>(height) * (1 + static_cast<size_t>(width) * 3));
        for (int y = 0; y < height; ++y) {
            raw.push_back(0); // filter type = None
            for (int x = 0; x < width; ++x) {
                size_t base = static_cast<size_t>(y * width + x) * 4;
                raw.push_back(pixels[base + 0]); // R
                raw.push_back(pixels[base + 1]); // G
                raw.push_back(pixels[base + 2]); // B
                // alpha ignored — color type 2
            }
        }
        auto compressed = zlib_stored(raw.data(), raw.size());
        write_chunk(out, "IDAT", compressed);
    }

    // IEND
    write_chunk(out, "IEND", {});

    return out;
}

// =============================================================================
// Rasterizer
// =============================================================================

static inline double clamp01(double v) {
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}

// Barycentric weights for point P w.r.t. triangle (A, B, C) in 2D screen-space.
// Returns false if degenerate.
static bool barycentric(double ax, double ay,
                         double bx, double by,
                         double cx, double cy,
                         double px, double py,
                         double& u, double& v, double& w) {
    double denom = (by - cy)*(ax - cx) + (cx - bx)*(ay - cy);
    if (std::fabs(denom) < 1e-10) return false;
    u = ((by - cy)*(px - cx) + (cx - bx)*(py - cy)) / denom;
    v = ((cy - ay)*(px - cx) + (ax - cx)*(py - cy)) / denom;
    w = 1.0 - u - v;
    return true;
}

std::vector<uint8_t> render_stl(
        const std::vector<float>& triangles,
        const Camera& cam,
        uint32_t bg_color,
        uint32_t mesh_color) {

    const int W = cam.width;
    const int H = cam.height;

    // Framebuffer: RGBA
    std::vector<uint8_t> fb(static_cast<size_t>(W * H) * 4, 0);
    // Z-buffer (stores -z in camera space, higher = closer)
    std::vector<double> zbuf(static_cast<size_t>(W * H),
                             -std::numeric_limits<double>::infinity());

    // Fill background
    uint8_t bg_r = (bg_color >> 16) & 0xFF;
    uint8_t bg_g = (bg_color >>  8) & 0xFF;
    uint8_t bg_b = (bg_color      ) & 0xFF;
    for (int i = 0; i < W * H; ++i) {
        fb[i*4+0] = bg_r;
        fb[i*4+1] = bg_g;
        fb[i*4+2] = bg_b;
        fb[i*4+3] = 255;
    }

    if (triangles.size() < 9) {
        return encode_png(fb, W, H);
    }

    // Camera setup
    Vec3 eye{cam.eye[0], cam.eye[1], cam.eye[2]};
    Vec3 tgt{cam.target[0], cam.target[1], cam.target[2]};
    Vec3 up {cam.up[0],  cam.up[1],  cam.up[2]};

    Mat4 view = look_at(eye, tgt, up);

    double fov_rad = cam.fov_deg * M_PI / 180.0;
    double aspect  = static_cast<double>(W) / H;
    // near/far: pick sensible values
    // First compute scene extent to set near/far
    double scene_r = 1.0;
    {
        // Approximate bounding sphere radius of mesh relative to centroid
        Vec3 cen{0,0,0};
        size_t n3 = triangles.size() / 3;
        for (size_t i = 0; i < n3; ++i)
            cen = cen + Vec3{triangles[i*3], triangles[i*3+1], triangles[i*3+2]};
        if (n3 > 0) cen = cen / static_cast<double>(n3);
        for (size_t i = 0; i < n3; ++i) {
            Vec3 p{triangles[i*3], triangles[i*3+1], triangles[i*3+2]};
            double d = (p - cen).len();
            if (d > scene_r) scene_r = d;
        }
    }
    double eye_dist = (eye - tgt).len();
    double near_z = std::max(0.001, eye_dist - scene_r * 2.0);
    double far_z  = eye_dist + scene_r * 4.0;

    Mat4 proj = perspective(fov_rad, aspect, near_z, far_z);
    Mat4 mvp  = proj * view;

    // Light direction = camera forward (from eye to target), in world space
    Vec3 light_dir = (tgt - eye).norm();
    // View direction = eye to fragment (for specular), simplified as light dir since same
    Vec3 view_dir  = light_dir;

    // Mesh base color
    double mr = ((mesh_color >> 16) & 0xFF) / 255.0;
    double mg = ((mesh_color >>  8) & 0xFF) / 255.0;
    double mb = ((mesh_color      ) & 0xFF) / 255.0;

    const double ka = 0.15;  // ambient
    const double kd = 0.70;  // diffuse
    const double ks = 0.15;  // specular
    const int    sh = 32;    // shininess

    size_t num_tris = triangles.size() / 9;

    for (size_t ti = 0; ti < num_tris; ++ti) {
        const float* base = triangles.data() + ti * 9;

        Vec3 v0{base[0], base[1], base[2]};
        Vec3 v1{base[3], base[4], base[5]};
        Vec3 v2{base[6], base[7], base[8]};

        // Face normal in world space
        Vec3 edge1 = v1 - v0;
        Vec3 edge2 = v2 - v0;
        Vec3 N = edge1.cross(edge2).norm();

        // Transform to NDC via MVP
        Vec3 c0 = mvp.mulPoint(v0);
        Vec3 c1 = mvp.mulPoint(v1);
        Vec3 c2 = mvp.mulPoint(v2);

        // Back-face culling in NDC (optional, improves speed)
        // Already handled by N dot view_dir check for shading

        // Convert NDC [-1,1] to screen pixels
        // NDC x: -1 = left, +1 = right
        // NDC y: -1 = bottom, +1 = top  → flip y for screen (top=0)
        auto ndc_to_px = [&](const Vec3& ndc, double& sx, double& sy) {
            sx = (ndc.x + 1.0) * 0.5 * W;
            sy = (1.0 - ndc.y) * 0.5 * H;
        };

        double sx0, sy0, sx1, sy1, sx2, sy2;
        ndc_to_px(c0, sx0, sy0);
        ndc_to_px(c1, sx1, sy1);
        ndc_to_px(c2, sx2, sy2);

        // Bounding box on screen
        int xmin = static_cast<int>(std::floor(std::min({sx0, sx1, sx2})));
        int xmax = static_cast<int>(std::ceil (std::max({sx0, sx1, sx2})));
        int ymin = static_cast<int>(std::floor(std::min({sy0, sy1, sy2})));
        int ymax = static_cast<int>(std::ceil (std::max({sy0, sy1, sy2})));

        xmin = std::max(xmin, 0);
        xmax = std::min(xmax, W - 1);
        ymin = std::max(ymin, 0);
        ymax = std::min(ymax, H - 1);

        if (xmin > xmax || ymin > ymax) continue;

        // Phong shading — compute per face (flat shading)
        double ndotl = N.dot(light_dir);
        double diffuse = kd * std::max(0.0, ndotl);

        // Specular: reflect light about normal
        Vec3 R = (N * (2.0 * ndotl) - light_dir).norm();
        double rdotv = R.dot(view_dir);
        double spec  = ks * std::pow(std::max(0.0, rdotv), sh);

        double shade = ka + diffuse + spec;
        // Clamp shade
        shade = std::min(shade, 1.2); // allow slight overexposure for highlights

        // Also try the other normal direction if facing away
        if (ndotl < 0.0) {
            // Flip and recompute for two-sided rendering
            Vec3 Nf = N * -1.0;
            double ndotl2 = Nf.dot(light_dir);
            double d2 = kd * std::max(0.0, ndotl2);
            Vec3 R2 = (Nf * (2.0 * ndotl2) - light_dir).norm();
            double s2 = ks * std::pow(std::max(0.0, R2.dot(view_dir)), sh);
            shade = ka + d2 + s2;
        }

        uint8_t pr = static_cast<uint8_t>(std::min(255.0, mr * shade * 255.0));
        uint8_t pg = static_cast<uint8_t>(std::min(255.0, mg * shade * 255.0));
        uint8_t pb = static_cast<uint8_t>(std::min(255.0, mb * shade * 255.0));

        for (int py = ymin; py <= ymax; ++py) {
            for (int px = xmin; px <= xmax; ++px) {
                double u, v, w;
                if (!barycentric(sx0, sy0, sx1, sy1, sx2, sy2,
                                 px + 0.5, py + 0.5, u, v, w))
                    continue;
                if (u < 0 || v < 0 || w < 0) continue;

                // Interpolate z (NDC z) for depth test
                double z = u * c0.z + v * c1.z + w * c2.z;

                // In NDC, z is in [-1,1]; closer = smaller z.
                // Store negated so that closer == larger value for max test.
                double depth = -z;

                size_t idx = static_cast<size_t>(py * W + px);
                if (depth <= zbuf[idx]) continue;
                zbuf[idx] = depth;

                fb[idx*4+0] = pr;
                fb[idx*4+1] = pg;
                fb[idx*4+2] = pb;
                fb[idx*4+3] = 255;
            }
        }
    }

    return encode_png(fb, W, H);
}

} // namespace render
