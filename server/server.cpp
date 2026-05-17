// =============================================================================
// CFD Mesh Server — HTTP/1.1 server integrating the CFD meshing pipeline.
// Compile from /Users/itkanmac2/cfd-mesh-cpp:
//   c++ -std=c++20 -O2 -I.. geometry/geometry_engine.cpp meshing/mesh_engine.cpp \
//       server/server.cpp -o cfd_server -lpthread
// =============================================================================

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <sys/stat.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstring>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "../geometry/geometry_engine.hpp"
#include "../geometry/mesh_sdf.hpp"
#include "../geometry/bifilar.hpp"
#include "../meshing/mesh_engine.hpp"
#include "../analysis/curvature.hpp"
#include "../slicer/slicer.hpp"
#include "../export/mesh_exporter.hpp"
#include "../geometry/primitives.hpp"
#include "../supabase/supabase_client.hpp"
#include "../shaders/shader_library.hpp"
#include "../shaders/node_compiler.hpp"
#include "../render/software_renderer.hpp"
#include "../cfd/lbm_sim.hpp"
#include "../cfd/flow3d.hpp"

using namespace cfd;

// =============================================================================
// In-memory stores
// =============================================================================

struct UploadRecord {
    std::string path;
    std::string filename;
    int         tri_count;
    std::vector<std::vector<uint32_t>>             face_groups;     // BFS co-planar triangle groups
    std::vector<std::vector<std::vector<Vec3d>>>   face_boundaries; // [group][loop][vertex] 3-D coords
};

struct MeshJob {
    std::string status;   // "pending","running","complete","error"
    std::string message;
    double      progress = 0.0;
    VolumeMesh  vol_mesh;
    PrismMesh   bl_mesh;
    double      mesh_size = 0.1;
    std::mutex  mu;
};

static std::map<std::string, UploadRecord>              uploads;
static std::map<std::string, std::shared_ptr<MeshJob>>  jobs;
static std::mutex                                        store_mu;

// Hot-reload: incremented by the file-watcher thread whenever webui files change.
static std::atomic<uint64_t> g_reload_gen{1};

// =============================================================================
// UUID / ID generation — 8 random hex characters
// =============================================================================

static std::string gen_id() {
    static std::random_device              rd;
    static std::mt19937                    rng(rd());
    static std::uniform_int_distribution<> dist(0, 15);
    static const char hex[] = "0123456789abcdef";
    std::string id;
    id.reserve(8);
    for (int i = 0; i < 8; ++i)
        id += hex[dist(rng)];
    return id;
}

// =============================================================================
// String utilities
// =============================================================================

static std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    size_t e = s.find_last_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    return s.substr(b, e - b + 1);
}

static bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

// Strip everything at and after '?' from a path segment.
static std::string strip_query(const std::string& s) {
    size_t q = s.find('?');
    return (q == std::string::npos) ? s : s.substr(0, q);
}

// Return everything after prefix, stripped of query string.
static std::string path_after(const std::string& path, const std::string& prefix) {
    if (!starts_with(path, prefix)) return "";
    return strip_query(path.substr(prefix.size()));
}

// =============================================================================
// JSON escaping
// =============================================================================

static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 4);
    for (unsigned char c : s) {
        if      (c == '"')  { out += "\\\""; }
        else if (c == '\\') { out += "\\\\"; }
        else if (c == '\n') { out += "\\n";  }
        else if (c == '\r') { out += "\\r";  }
        else if (c == '\t') { out += "\\t";  }
        else if (c < 0x20)  { out += "\\u00";
                               out += "0123456789abcdef"[c >> 4];
                               out += "0123456789abcdef"[c & 0xf]; }
        else                { out += static_cast<char>(c); }
    }
    return out;
}

// =============================================================================
// HTTP response builders
// =============================================================================

static std::string http_response(int code, const std::string& reason,
                                 const std::string& content_type,
                                 const std::string& body) {
    std::ostringstream os;
    os << "HTTP/1.1 " << code << " " << reason << "\r\n"
       << "Content-Type: " << content_type << "\r\n"
       << "Content-Length: " << body.size() << "\r\n"
       << "Connection: close\r\n"
       << "Access-Control-Allow-Origin: *\r\n"
       << "Cache-Control: no-store\r\n"
       << "\r\n"
       << body;
    return os.str();
}

static std::string ok_json(const std::string& body) {
    return http_response(200, "OK", "application/json", body);
}

static std::string accepted_json(const std::string& body) {
    return http_response(202, "Accepted", "application/json", body);
}

static std::string error_json(int code, const std::string& msg) {
    const char* reason = (code == 404) ? "Not Found"   :
                         (code == 400) ? "Bad Request"  :
                         (code == 405) ? "Method Not Allowed" :
                                         "Internal Server Error";
    std::string body = "{\"error\":\"" + json_escape(msg) + "\"}";
    return http_response(code, reason, "application/json", body);
}

static std::string mime_for(const std::string& path) {
    size_t dot = path.rfind('.');
    if (dot == std::string::npos) return "application/octet-stream";
    std::string ext = path.substr(dot);
    if (ext == ".html") return "text/html; charset=utf-8";
    if (ext == ".css")  return "text/css";
    if (ext == ".js")   return "application/javascript";
    if (ext == ".json") return "application/json";
    if (ext == ".svg")  return "image/svg+xml";
    if (ext == ".png")  return "image/png";
    if (ext == ".jpg" || ext == ".jpeg") return "image/jpeg";
    if (ext == ".ico")  return "image/x-icon";
    if (ext == ".woff") return "font/woff";
    if (ext == ".woff2") return "font/woff2";
    if (ext == ".ttf")  return "font/ttf";
    if (ext == ".eot")  return "application/vnd.ms-fontobject";
    if (ext == ".stl")  return "model/stl";
    return "application/octet-stream";
}

static std::string serve_file(const std::string& path, const std::string& mime) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open())
        return error_json(404, "File not found: " + path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return http_response(200, "OK", mime, ss.str());
}

// =============================================================================
// HTTP request parser
// =============================================================================

struct HttpRequest {
    std::string                       method;
    std::string                       path;
    std::string                       version;
    std::map<std::string,std::string> headers;
    std::string                       body;
    bool                              valid = false;
};

// Read from socket until we have a complete HTTP request (headers + full body).
// Handles large uploads up to 50 MB without re-parsing headers on each chunk.
static std::string recv_all(int fd) {
    std::string buf;
    buf.reserve(65536);
    char tmp[65536];  // Larger chunk size: fewer iterations for big uploads

    // Phase 1: read until we have the complete header block (\r\n\r\n).
    size_t hdr_end = std::string::npos;
    while (hdr_end == std::string::npos) {
        ssize_t n = recv(fd, tmp, sizeof(tmp), 0);
        if (n <= 0) return buf;
        buf.append(tmp, static_cast<size_t>(n));
        hdr_end = buf.find("\r\n\r\n");
    }

    // Phase 2: parse Content-Length exactly once from the already-received headers.
    std::string lower_hdr = to_lower(buf.substr(0, hdr_end));
    const std::string cl_key = "content-length:";
    size_t cl_pos = lower_hdr.find(cl_key);
    if (cl_pos == std::string::npos) return buf; // no body expected

    size_t val_start = cl_pos + cl_key.size();
    size_t val_end   = lower_hdr.find('\n', val_start);
    long content_length = 0;
    try {
        std::string cl_str = (val_end == std::string::npos)
                             ? lower_hdr.substr(val_start)
                             : lower_hdr.substr(val_start, val_end - val_start);
        content_length = std::stol(cl_str);
    } catch (...) { return buf; }

    if (content_length <= 0) return buf;

    // Pre-allocate for the full request to avoid repeated reallocations on
    // large uploads (up to 50 MB).
    size_t body_start = hdr_end + 4;
    size_t total_expected = body_start + static_cast<size_t>(content_length);
    if (total_expected > buf.capacity())
        buf.reserve(total_expected);

    // Phase 3: keep reading until we have all body bytes.
    while (static_cast<long>(buf.size() - body_start) < content_length) {
        ssize_t n = recv(fd, tmp, sizeof(tmp), 0);
        if (n <= 0) break;
        buf.append(tmp, static_cast<size_t>(n));
    }
    return buf;
}

static HttpRequest parse_request(const std::string& raw) {
    HttpRequest req;
    if (raw.empty()) return req;

    size_t hdr_end = raw.find("\r\n\r\n");
    if (hdr_end == std::string::npos) return req;

    std::string hdr_section = raw.substr(0, hdr_end);
    req.body = raw.substr(hdr_end + 4);

    // Request line
    size_t line_end = hdr_section.find("\r\n");
    std::string req_line = (line_end == std::string::npos)
                           ? hdr_section : hdr_section.substr(0, line_end);
    {
        std::istringstream rls(req_line);
        rls >> req.method >> req.path >> req.version;
    }
    if (req.method.empty() || req.path.empty()) return req;

    // Headers
    size_t pos = (line_end == std::string::npos) ? hdr_section.size() : line_end + 2;
    while (pos < hdr_section.size()) {
        size_t next = hdr_section.find("\r\n", pos);
        size_t len  = (next == std::string::npos) ? std::string::npos : next - pos;
        std::string line = hdr_section.substr(pos, len);
        size_t colon = line.find(':');
        if (colon != std::string::npos) {
            std::string key = to_lower(trim(line.substr(0, colon)));
            std::string val = trim(line.substr(colon + 1));
            req.headers[key] = val;
        }
        if (next == std::string::npos) break;
        pos = next + 2;
    }

    // Enforce Content-Length on body
    auto cl_it = req.headers.find("content-length");
    if (cl_it != req.headers.end()) {
        try {
            size_t cl = std::stoul(cl_it->second);
            if (req.body.size() > cl) req.body.resize(cl);
        } catch (...) {}
    }

    req.valid = true;
    return req;
}

// =============================================================================
// Multipart/form-data parser
// =============================================================================

struct MultipartPart {
    std::map<std::string, std::string> headers;
    std::string                        data;
    std::string                        name;
    std::string                        filename;
};

static std::string extract_boundary(const std::string& ct) {
    const std::string key = "boundary=";
    size_t pos = ct.find(key);
    if (pos == std::string::npos) return "";
    std::string bnd = ct.substr(pos + key.size());
    // Strip quotes
    if (!bnd.empty() && bnd.front() == '"') {
        bnd = bnd.substr(1);
        size_t q = bnd.find('"');
        if (q != std::string::npos) bnd.resize(q);
    }
    // Stop at semicolon or whitespace
    size_t stop = bnd.find_first_of("; \t\r\n");
    if (stop != std::string::npos) bnd.resize(stop);
    return trim(bnd);
}

static std::string extract_param(const std::string& hdr, const std::string& param_name) {
    std::string k = param_name + "=";
    size_t pos = hdr.find(k);
    if (pos == std::string::npos) return "";
    pos += k.size();
    if (pos >= hdr.size()) return "";
    if (hdr[pos] == '"') {
        ++pos;
        size_t end = hdr.find('"', pos);
        if (end == std::string::npos) return hdr.substr(pos);
        return hdr.substr(pos, end - pos);
    }
    size_t end = hdr.find_first_of(";\r\n ", pos);
    if (end == std::string::npos) return hdr.substr(pos);
    return hdr.substr(pos, end - pos);
}

static std::vector<MultipartPart> parse_multipart(const std::string& body,
                                                   const std::string& boundary) {
    std::vector<MultipartPart> parts;
    const std::string delim = "--" + boundary;

    size_t pos = body.find(delim);
    if (pos == std::string::npos) return parts;

    while (true) {
        pos += delim.size();

        // Check for terminator "--"
        if (pos + 1 < body.size() && body[pos] == '-' && body[pos+1] == '-')
            break;
        // Skip CRLF after boundary line
        if (pos < body.size() && body[pos] == '\r') ++pos;
        if (pos < body.size() && body[pos] == '\n') ++pos;

        MultipartPart part;

        // Parse part headers until blank line
        while (true) {
            size_t eol = body.find("\r\n", pos);
            if (eol == std::string::npos) return parts; // malformed
            std::string line = body.substr(pos, eol - pos);
            pos = eol + 2;
            if (line.empty()) break; // blank line ends headers
            size_t colon = line.find(':');
            if (colon != std::string::npos) {
                std::string hk = to_lower(trim(line.substr(0, colon)));
                std::string hv = trim(line.substr(colon + 1));
                part.headers[hk] = hv;
            }
        }

        // Extract name / filename from Content-Disposition
        {
            auto it = part.headers.find("content-disposition");
            if (it != part.headers.end()) {
                part.name     = extract_param(it->second, "name");
                part.filename = extract_param(it->second, "filename");
            }
        }

        // Data ends at the next boundary (preceded by CRLF)
        {
            std::string next_delim = "\r\n" + delim;
            size_t data_end = body.find(next_delim, pos);
            if (data_end == std::string::npos) {
                // Last part — try without leading CRLF
                data_end = body.find(delim, pos);
                if (data_end == std::string::npos) data_end = body.size();
                part.data = body.substr(pos, data_end - pos);
                parts.push_back(std::move(part));
                return parts;
            }
            part.data = body.substr(pos, data_end - pos);
            parts.push_back(std::move(part));
            pos = data_end + 2; // advance past the CRLF before next --boundary
        }
    }
    return parts;
}

// =============================================================================
// Minimal JSON value extraction (no external library)
// =============================================================================

static std::string json_string_value(const std::string& json, const std::string& key) {
    std::string search = "\"" + key + "\"";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return "";
    pos += search.size();
    // Skip whitespace and colon
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == ':' || json[pos] == '\t'))
        ++pos;
    if (pos >= json.size() || json[pos] != '"') return "";
    ++pos;
    std::string val;
    while (pos < json.size() && json[pos] != '"') {
        if (json[pos] == '\\' && pos + 1 < json.size()) {
            char esc = json[++pos];
            if      (esc == '"')  { val += '"';  }
            else if (esc == '\\') { val += '\\'; }
            else if (esc == 'n')  { val += '\n'; }
            else if (esc == 'r')  { val += '\r'; }
            else if (esc == 't')  { val += '\t'; }
            else                  { val += esc;  }
            ++pos;
        } else {
            val += json[pos++];
        }
    }
    return val;
}

static double json_double_value(const std::string& json, const std::string& key,
                                double def = 0.0) {
    std::string search = "\"" + key + "\"";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return def;
    pos += search.size();
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == ':' || json[pos] == '\t'))
        ++pos;
    if (pos >= json.size()) return def;
    try {
        size_t consumed = 0;
        double v = std::stod(json.substr(pos), &consumed);
        return (consumed > 0) ? v : def;
    } catch (...) { return def; }
}

static int json_int_value(const std::string& json, const std::string& key, int def = 0) {
    return static_cast<int>(json_double_value(json, key, static_cast<double>(def)));
}

// =============================================================================
// Mesh data helpers
// =============================================================================

// A boundary face of a tet is identified by its 3 sorted vertex indices.
using Face3 = std::array<uint32_t, 3>;
using Edge2 = std::array<uint32_t, 2>;

static Face3 make_face(uint32_t a, uint32_t b, uint32_t c) {
    Face3 f = {a, b, c};
    std::sort(f.begin(), f.end());
    return f;
}

// All 4 faces of a tetrahedron.
static std::array<Face3, 4> tet_faces(const TetElement& t) {
    const auto& n = t.nodes;
    return {{
        make_face(n[0], n[1], n[2]),
        make_face(n[0], n[1], n[3]),
        make_face(n[0], n[2], n[3]),
        make_face(n[1], n[2], n[3])
    }};
}

static Edge2 make_edge(uint32_t a, uint32_t b) {
    return (a < b) ? Edge2{a, b} : Edge2{b, a};
}

// =============================================================================
// Fallback STL triangle counter
//
// Used when GeometryEngine::importCAD throws. Does not depend on the meshing
// pipeline at all — just reads the raw bytes to count faces.
//
// ASCII STL:  count occurrences of the keyword "facet normal"
// Binary STL: read the 4-byte triangle count at offset 80 and validate it
//             against the file size (80 header + 4 count + n*50 bytes).
// =============================================================================

static int count_stl_triangles(const std::string& data) {
    if (data.size() < 84) {
        // Could be a very small ASCII STL; fall through to ASCII scan.
    }

    // Decide ASCII vs binary the same way GeometryEngine does:
    // if first 5 bytes == "solid" it *might* be ASCII, but verify via size.
    bool try_binary = true;
    if (data.size() >= 5 && data.compare(0, 5, "solid") == 0) {
        // Could be ASCII. Check the binary consistency.
        if (data.size() >= 84) {
            uint32_t bin_count = 0;
            std::memcpy(&bin_count, data.data() + 80, 4);
            size_t expected = 84u + static_cast<size_t>(bin_count) * 50u;
            try_binary = (data.size() == expected) && (bin_count > 0);
        } else {
            try_binary = false;
        }
    }

    if (try_binary && data.size() >= 84) {
        uint32_t bin_count = 0;
        std::memcpy(&bin_count, data.data() + 80, 4);
        size_t expected = 84u + static_cast<size_t>(bin_count) * 50u;
        if (data.size() == expected && bin_count > 0) {
            return static_cast<int>(bin_count);
        }
    }

    // ASCII fallback: count "facet normal" substrings
    int count = 0;
    const std::string needle = "facet normal";
    size_t pos = 0;
    while ((pos = data.find(needle, pos)) != std::string::npos) {
        ++count;
        pos += needle.size();
    }
    return count;
}

// =============================================================================
// API handlers
// =============================================================================

static std::string handle_health() {
    return ok_json("{\"status\":\"ok\"}");
}

// POST /api/upload  — multipart/form-data with a "file" part
static std::string handle_upload(const HttpRequest& req) {
    auto ct_it = req.headers.find("content-type");
    if (ct_it == req.headers.end())
        return error_json(400, "Missing Content-Type header");

    std::string boundary = extract_boundary(ct_it->second);
    if (boundary.empty())
        return error_json(400, "Cannot parse multipart boundary from Content-Type");

    auto parts = parse_multipart(req.body, boundary);

    MultipartPart* file_part = nullptr;
    for (auto& p : parts) {
        if (p.name == "file") { file_part = &p; break; }
    }
    if (!file_part)
        return error_json(400, "No field named 'file' in multipart body");
    if (file_part->data.empty())
        return error_json(400, "Uploaded file is empty");

    // Write to /tmp/<id>.stl
    std::string id       = gen_id();
    std::string tmp_path = "/tmp/" + id + ".stl";
    {
        std::ofstream out(tmp_path, std::ios::binary);
        if (!out.is_open())
            return error_json(500, "Cannot open " + tmp_path + " for writing");
        out.write(file_part->data.data(), static_cast<std::streamsize>(file_part->data.size()));
        if (!out)
            return error_json(500, "Write to " + tmp_path + " failed");
    }

    // Import geometry, merge vertices, and compute BFS flat-face groups.
    // Fall back to a lightweight triangle counter if the geometry engine fails.
    int tri_count = 0;
    std::vector<std::vector<uint32_t>>           face_groups;
    std::vector<std::vector<std::vector<Vec3d>>> face_boundaries;
    try {
        GeometryEngine geo([](const std::string&, double) {});
        SurfaceManifold sm = geo.importCAD(tmp_path);
        geo.mergeVertices(sm);
        face_groups     = geo.identifyFlatFaces(sm);
        auto loopIdxs   = geo.extractAllBoundaryLoops(sm, face_groups, /*simplify=*/true);
        tri_count       = static_cast<int>(sm.faces.size());

        // Convert vertex indices → 3-D coordinates for the API response
        face_boundaries.reserve(loopIdxs.size());
        for (const auto& loops : loopIdxs) {
            std::vector<std::vector<Vec3d>> loopCoords;
            loopCoords.reserve(loops.size());
            for (const auto& loop : loops) {
                std::vector<Vec3d> coords;
                coords.reserve(loop.size());
                for (uint32_t vi : loop) {
                    coords.push_back(sm.vertices[vi].pos);
                }
                loopCoords.push_back(std::move(coords));
            }
            face_boundaries.push_back(std::move(loopCoords));
        }
    } catch (const std::exception&) {
        tri_count = count_stl_triangles(file_part->data);
    }

    std::string filename = file_part->filename.empty() ? "upload.stl" : file_part->filename;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        uploads[id] = UploadRecord{
            tmp_path, filename, tri_count,
            std::move(face_groups), std::move(face_boundaries)
        };
    }

    std::ostringstream resp;
    resp << "{\"file_id\":\"" << id
         << "\",\"triangle_count\":" << tri_count
         << ",\"filename\":\"" << json_escape(filename) << "\"}";
    return ok_json(resp.str());
}

// POST /api/generate-mesh
static std::string handle_generate_mesh(const HttpRequest& req) {
    const std::string& body = req.body;

    std::string file_id = json_string_value(body, "file_id");
    double mesh_size    = json_double_value(body, "mesh_size",       0.1);
    int    bl_layers    = json_int_value   (body, "boundary_layers", 5);
    double bl_ratio     = json_double_value(body, "bl_ratio",        1.2);

    if (file_id.empty())
        return error_json(400, "Missing required field: file_id");

    UploadRecord rec;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(file_id);
        if (it == uploads.end())
            return error_json(404, "Unknown file_id: " + file_id);
        rec = it->second;
    }

    std::string mesh_id = gen_id();
    auto job = std::make_shared<MeshJob>();
    job->status    = "pending";
    job->message   = "Queued";
    job->mesh_size = mesh_size;

    {
        std::lock_guard<std::mutex> lk(store_mu);
        jobs[mesh_id] = job;
    }

    // Launch the pipeline in a background thread
    std::thread([job, rec, mesh_size, bl_layers, bl_ratio]() mutable {

        auto set_progress = [&](const std::string& msg, double pct) {
            std::lock_guard<std::mutex> lk(job->mu);
            job->message  = msg;
            job->progress = pct;
        };

        try {
            {
                std::lock_guard<std::mutex> lk(job->mu);
                job->status   = "running";
                job->message  = "Importing geometry";
                job->progress = 0.0;
            }

            PipelineConfig cfg;
            cfg.input_file           = rec.path;
            cfg.reynolds_number      = 1e6;
            cfg.target_yplus         = 1.0;
            cfg.kinematic_viscosity  = 1.5e-5;
            cfg.freestream_velocity  = 10.0;
            cfg.bl_layers            = bl_layers;
            cfg.bl_growth_ratio      = bl_ratio;
            cfg.core_mesh_size       = mesh_size;
            cfg.output_format        = "openfoam";

            // --- Geometry stage (0 → 30%) ---
            GeometryEngine geo([&](const std::string& stage, double pct) {
                set_progress("Geometry: " + stage, pct * 0.30);
            });
            SurfaceManifold sm = geo.importCAD(rec.path);
            set_progress("Repairing surface", 0.25);
            (void)geo.repairManifold(sm);
            geo.identifyBoundaries(sm, {});   // tag all unlabelled faces as Wall
            set_progress("Surface ready", 0.30);

            // --- Meshing stage (30% → 100%) ---
            MeshEngine mesher([&](const std::string& stage, double pct) {
                set_progress("Meshing: " + stage, 0.30 + pct * 0.70);
            });

            BLParameters bl_params = MeshEngine::computeBLParameters(cfg);
            set_progress("Inflating prisms", 0.35);

            PrismMesh  bl  = mesher.inflatePrisms(sm, bl_params);
            set_progress("Generating tet core", 0.65);

            VolumeMesh vol = mesher.generateTetCore(bl, mesh_size);
            set_progress("Checking quality", 0.95);

            (void)mesher.checkQuality(vol);

            {
                std::lock_guard<std::mutex> lk(job->mu);
                job->vol_mesh = std::move(vol);
                job->bl_mesh  = std::move(bl);
                job->status   = "complete";
                job->message  = "Mesh generation complete";
                job->progress = 1.0;
            }

        } catch (const std::exception& ex) {
            std::lock_guard<std::mutex> lk(job->mu);
            job->status   = "error";
            job->message  = ex.what();
        } catch (...) {
            std::lock_guard<std::mutex> lk(job->mu);
            job->status  = "error";
            job->message = "Unknown error during mesh generation";
        }

    }).detach();

    std::string resp = "{\"mesh_id\":\"" + mesh_id + "\"}";
    return accepted_json(resp);
}

// GET /api/mesh/<id>/status
static std::string handle_mesh_status(const std::string& mesh_id) {
    std::shared_ptr<MeshJob> job;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = jobs.find(mesh_id);
        if (it == jobs.end())
            return error_json(404, "Unknown mesh_id: " + mesh_id);
        job = it->second;
    }

    std::string jstatus, jmessage;
    double jprogress;
    {
        std::lock_guard<std::mutex> lk(job->mu);
        jstatus   = job->status;
        jmessage  = job->message;
        jprogress = job->progress;
    }

    std::ostringstream resp;
    resp << "{\"status\":\"" << jstatus
         << "\",\"progress\":" << jprogress
         << ",\"message\":\"" << json_escape(jmessage) << "\"}";
    return ok_json(resp.str());
}

// GET /api/mesh/<id>/data
static std::string handle_mesh_data(const std::string& mesh_id) {
    std::shared_ptr<MeshJob> job;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = jobs.find(mesh_id);
        if (it == jobs.end())
            return error_json(404, "Unknown mesh_id: " + mesh_id);
        job = it->second;
    }

    VolumeMesh vol;
    PrismMesh  bl;
    {
        std::lock_guard<std::mutex> lk(job->mu);
        if (job->status != "complete")
            return error_json(400, "Mesh not ready (status: " + job->status + ")");
        vol = job->vol_mesh;
        bl  = job->bl_mesh;
    }

    std::ostringstream out;
    out << "{";

    // --- nodes ---
    out << "\"nodes\":[";
    for (size_t i = 0; i < vol.vertices.size(); ++i) {
        if (i) out << ",";
        const auto& p = vol.vertices[i].pos;
        out << "[" << p[0] << "," << p[1] << "," << p[2] << "]";
    }
    out << "]";

    // --- surface_triangles ---
    // Faces that appear in exactly one tet are on the surface boundary.
    int total_surface_tris = 0;
    {
        std::map<Face3, int> face_count;
        for (const auto& tet : vol.tets)
            for (const auto& f : tet_faces(tet))
                face_count[f]++;

        out << ",\"surface_triangles\":[";
        bool first = true;
        int  emitted = 0;
        for (const auto& kv : face_count) {
            if (kv.second != 1) continue;
            ++total_surface_tris;
            if (emitted < 20000) {
                const Face3& f = kv.first;
                if (!first) out << ",";
                out << "[" << f[0] << "," << f[1] << "," << f[2] << "]";
                first = false;
                ++emitted;
            }
        }
        out << "]";
    }

    // --- volume_edges (unique edges from tets, limit 50000) ---
    {
        std::set<Edge2> es;
        for (const auto& tet : vol.tets) {
            const auto& n = tet.nodes;
            for (int i = 0; i < 4 && static_cast<int>(es.size()) < 50000; ++i)
                for (int j = i+1; j < 4 && static_cast<int>(es.size()) < 50000; ++j)
                    es.insert(make_edge(n[i], n[j]));
            if (static_cast<int>(es.size()) >= 50000) break;
        }
        out << ",\"volume_edges\":[";
        bool first = true;
        for (const auto& e : es) {
            if (!first) out << ",";
            out << "[" << e[0] << "," << e[1] << "]";
            first = false;
        }
        out << "]";
    }

    // --- bl_edges (edges from prism elements in vol.prisms, limit 20000) ---
    {
        std::set<Edge2> be;
        for (const auto& prism : vol.prisms) {
            const auto& n = prism.nodes;
            // Bottom triangle: 0-1-2
            be.insert(make_edge(n[0], n[1]));
            be.insert(make_edge(n[1], n[2]));
            be.insert(make_edge(n[0], n[2]));
            // Top triangle: 3-4-5
            be.insert(make_edge(n[3], n[4]));
            be.insert(make_edge(n[4], n[5]));
            be.insert(make_edge(n[3], n[5]));
            // Lateral edges: 0-3, 1-4, 2-5
            be.insert(make_edge(n[0], n[3]));
            be.insert(make_edge(n[1], n[4]));
            be.insert(make_edge(n[2], n[5]));
            if (static_cast<int>(be.size()) >= 20000) break;
        }
        out << ",\"bl_edges\":[";
        bool first = true;
        for (const auto& e : be) {
            if (!first) out << ",";
            out << "[" << e[0] << "," << e[1] << "]";
            first = false;
        }
        out << "]";
    }

    // --- bounds ---
    {
        double xmn =  1e18, ymn =  1e18, zmn =  1e18;
        double xmx = -1e18, ymx = -1e18, zmx = -1e18;
        for (const auto& vtx : vol.vertices) {
            xmn = std::min(xmn, vtx.pos[0]); xmx = std::max(xmx, vtx.pos[0]);
            ymn = std::min(ymn, vtx.pos[1]); ymx = std::max(ymx, vtx.pos[1]);
            zmn = std::min(zmn, vtx.pos[2]); zmx = std::max(zmx, vtx.pos[2]);
        }
        if (vol.vertices.empty()) { xmn=xmx=ymn=ymx=zmn=zmx=0.0; }
        double dx = xmx-xmn, dy = ymx-ymn, dz = zmx-zmn;
        double diag = std::sqrt(dx*dx + dy*dy + dz*dz);
        out << ",\"bounds\":{\"min\":["  << xmn << "," << ymn << "," << zmn
            << "],\"max\":["             << xmx << "," << ymx << "," << zmx
            << "],\"diagonal\":"         << diag << "}";
    }

    // --- stats ---
    out << ",\"stats\":{"
        << "\"num_nodes\":"       << vol.vertices.size()
        << ",\"num_tets\":"       << vol.tets.size()
        << ",\"num_surface_tris\":" << std::min(total_surface_tris, 20000)
        << ",\"num_prisms\":"     << vol.prisms.size()
        << "}";

    out << "}";
    return ok_json(out.str());
}

// GET /api/mesh/<id>/stl
// Returns surface preview from boundary layer prism outer faces (first 200 prisms).
static std::string handle_mesh_stl(const std::string& mesh_id) {
    std::shared_ptr<MeshJob> job;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = jobs.find(mesh_id);
        if (it == jobs.end())
            return error_json(404, "Unknown mesh_id: " + mesh_id);
        job = it->second;
    }

    PrismMesh bl;
    {
        std::lock_guard<std::mutex> lk(job->mu);
        if (job->status != "complete")
            return error_json(400, "Mesh not ready (status: " + job->status + ")");
        bl = job->bl_mesh;
    }

    int prism_limit = std::min(static_cast<int>(bl.prisms.size()), 200);

    // Collect the subset of nodes referenced by top faces (nodes[3..5]) of
    // the first prism_limit prisms, remapping to a compact index.
    std::map<uint32_t, uint32_t> node_remap;
    std::vector<uint32_t>        node_list;

    for (int pi = 0; pi < prism_limit; ++pi) {
        const auto& prism = bl.prisms[pi];
        for (int k = 3; k < 6; ++k) {
            uint32_t orig = prism.nodes[k];
            if (node_remap.find(orig) == node_remap.end()) {
                node_remap[orig] = static_cast<uint32_t>(node_list.size());
                node_list.push_back(orig);
            }
        }
    }

    std::ostringstream out;
    out << "{\"nodes\":[";
    for (size_t i = 0; i < node_list.size(); ++i) {
        if (i) out << ",";
        uint32_t idx = node_list[i];
        if (idx < bl.vertices.size()) {
            const auto& p = bl.vertices[idx].pos;
            out << "[" << p[0] << "," << p[1] << "," << p[2] << "]";
        } else {
            out << "[0,0,0]";
        }
    }
    out << "],\"triangles\":[";
    for (int pi = 0; pi < prism_limit; ++pi) {
        if (pi) out << ",";
        const auto& prism = bl.prisms[pi];
        uint32_t a = node_remap.at(prism.nodes[3]);
        uint32_t b = node_remap.at(prism.nodes[4]);
        uint32_t c = node_remap.at(prism.nodes[5]);
        out << "[" << a << "," << b << "," << c << "]";
    }
    out << "]}";

    return ok_json(out.str());
}

// GET /api/upload/<id>/stl
// Serves the raw STL file for a given upload_id.
static std::string handle_upload_stl(const std::string& upload_id) {
    std::string path, filename;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(upload_id);
        if (it == uploads.end())
            return error_json(404, "Unknown upload_id: " + upload_id);
        path     = it->second.path;
        filename = it->second.filename;
    }
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open())
        return error_json(404, "File not found on disk: " + path);
    std::string body((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());
    std::ostringstream hdr;
    hdr << "HTTP/1.1 200 OK\r\n"
        << "Content-Type: model/stl\r\n"
        << "Content-Length: " << body.size() << "\r\n"
        << "Connection: close\r\n"
        << "Access-Control-Allow-Origin: *\r\n"
        << "\r\n"
        << body;
    return hdr.str();
}

// GET /api/upload/<id>/faces
// Returns BFS flat-face groups computed at upload time.
// Response: { triangle_count, group_count, groups: [[tri_idx,...], ...] }
static std::string handle_upload_faces(const std::string& upload_id) {
    int tri_count = 0;
    std::vector<std::vector<uint32_t>>           groups;
    std::vector<std::vector<std::vector<Vec3d>>> boundaries;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(upload_id);
        if (it == uploads.end())
            return error_json(404, "Unknown upload_id: " + upload_id);
        tri_count  = it->second.tri_count;
        groups     = it->second.face_groups;
        boundaries = it->second.face_boundaries;
    }

    std::ostringstream out;
    out << std::fixed;
    out.precision(7);

    out << "{\"triangle_count\":" << tri_count
        << ",\"group_count\":"    << groups.size()
        << ",\"groups\":[";

    for (size_t gi = 0; gi < groups.size(); ++gi) {
        if (gi) out << ",";
        out << "[";
        const auto& grp = groups[gi];
        for (size_t j = 0; j < grp.size(); ++j) {
            if (j) out << ",";
            out << grp[j];
        }
        out << "]";
    }

    // boundaries[gi][loop_idx] = list of [x,y,z] coords
    out << "],\"boundaries\":[";
    for (size_t gi = 0; gi < boundaries.size(); ++gi) {
        if (gi) out << ",";
        out << "[";
        const auto& loops = boundaries[gi];
        for (size_t li = 0; li < loops.size(); ++li) {
            if (li) out << ",";
            out << "[";
            const auto& loop = loops[li];
            for (size_t vi = 0; vi < loop.size(); ++vi) {
                if (vi) out << ",";
                out << "[" << loop[vi][0] << ","
                           << loop[vi][1] << ","
                           << loop[vi][2] << "]";
            }
            out << "]";
        }
        out << "]";
    }
    out << "]}";

    return ok_json(out.str());
}

// POST /api/heal
static std::string handle_heal(const HttpRequest& req) {
    std::string file_id = json_string_value(req.body, "file_id");
    if (file_id.empty())
        return error_json(400, "Missing required field: file_id");

    std::string path, orig_name;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(file_id);
        if (it == uploads.end())
            return error_json(404, "Unknown file_id: " + file_id);
        path      = it->second.path;
        orig_name = it->second.filename;
    }

    try {
        GeometryEngine geo([](const std::string&, double) {});
        SurfaceManifold sm = geo.importCAD(path);
        bool watertight = geo.repairManifold(sm);
        auto report     = geo.analyzeManifold(sm);

        std::string new_id   = gen_id();
        std::string out_path = "/tmp/" + new_id + "_healed.stl";
        cfd::makeExporter("stl_binary")->write(sm, out_path);

        size_t dot  = orig_name.rfind('.');
        std::string stem = (dot == std::string::npos) ? orig_name : orig_name.substr(0, dot);
        {
            std::lock_guard<std::mutex> lk(store_mu);
            uploads[new_id] = UploadRecord{
                out_path, stem + "_healed.stl",
                static_cast<int>(sm.faces.size()),
                {}, {}
            };
        }

        std::ostringstream r;
        r << "{\"healed_file_id\":\"" << new_id << "\""
          << ",\"vertex_count\":"       << sm.vertices.size()
          << ",\"face_count\":"         << sm.faces.size()
          << ",\"is_watertight\":"      << (watertight ? "true" : "false")
          << ",\"boundary_edges\":"     << report.boundary_edges
          << ",\"non_manifold_edges\":" << report.non_manifold_edges
          << ",\"isolated_vertices\":"  << report.isolated_vertices
          << "}";
        return ok_json(r.str());
    } catch (const std::exception& e) {
        return error_json(500, std::string("Heal failed: ") + e.what());
    }
}

// POST /api/slice
static std::string handle_slice(const HttpRequest& req) {
    std::string file_id = json_string_value(req.body, "file_id");
    if (file_id.empty())
        return error_json(400, "Missing required field: file_id");
    double z = json_double_value(req.body, "z", 0.0);

    std::string path;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(file_id);
        if (it == uploads.end())
            return error_json(404, "Unknown file_id: " + file_id);
        path = it->second.path;
    }

    try {
        GeometryEngine geo([](const std::string&, double) {});
        SurfaceManifold sm = geo.importCAD(path);
        (void)geo.repairManifold(sm);

        double zmin =  1e18, zmax = -1e18;
        for (const auto& v : sm.vertices) {
            if (v.pos[2] < zmin) zmin = v.pos[2];
            if (v.pos[2] > zmax) zmax = v.pos[2];
        }

        cfd::Slicer slicer;
        auto result = slicer.sliceAtZ(sm, z);

        std::ostringstream r;
        r << std::fixed;
        r.precision(8);
        r << "{\"z_level\":" << z
          << ",\"z_min\":"   << zmin
          << ",\"z_max\":"   << zmax
          << ",\"contour_count\":" << result.contours.size()
          << ",\"contours\":[";

        bool fc = true;
        for (const auto& c : result.contours) {
            if (!fc) r << ",";
            r << "[";
            bool fp = true;
            for (const auto& p : c.points) {
                if (!fp) r << ",";
                r << "[" << p[0] << "," << p[1] << "]";
                fp = false;
            }
            r << "]";
            fc = false;
        }
        r << "]}";
        return ok_json(r.str());
    } catch (const std::exception& e) {
        return error_json(500, std::string("Slice failed: ") + e.what());
    }
}

// GET /api/curvature/<upload_id>
static std::string handle_curvature(const std::string& upload_id) {
    std::string path;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(upload_id);
        if (it == uploads.end())
            return error_json(404, "Unknown upload_id: " + upload_id);
        path = it->second.path;
    }

    try {
        GeometryEngine geo([](const std::string&, double) {});
        SurfaceManifold sm = geo.importCAD(path);
        (void)geo.repairManifold(sm);

        auto result        = cfd::CurvatureAnalyzer().compute(sm);
        const auto& H      = result.mean_curvature;
        const auto& jc     = result.jet_colors;
        const uint32_t nv  = static_cast<uint32_t>(sm.vertices.size());

        double hmin = *std::min_element(H.begin(), H.end());
        double hmax = *std::max_element(H.begin(), H.end());

        std::ostringstream r;
        r << std::fixed;
        r.precision(7);
        r << "{\"vertex_count\":" << nv
          << ",\"min_h\":"        << hmin
          << ",\"max_h\":"        << hmax
          << ",\"positions\":[";
        for (uint32_t i = 0; i < nv; ++i) {
            if (i) r << ",";
            const auto& p = sm.vertices[i].pos;
            r << p[0] << "," << p[1] << "," << p[2];
        }
        r << "],\"colors\":[";
        r.precision(5);
        for (uint32_t i = 0; i < nv; ++i) {
            if (i) r << ",";
            r << jc[i][0] / 255.0f << ","
              << jc[i][1] / 255.0f << ","
              << jc[i][2] / 255.0f;
        }
        r << "]}";
        return ok_json(r.str());
    } catch (const std::exception& e) {
        return error_json(500, std::string("Curvature failed: ") + e.what());
    }
}

// GET /api/export/<upload_id>/<format>
static std::string handle_export(const std::string& upload_id, const std::string& format) {
    std::string path, filename;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(upload_id);
        if (it == uploads.end())
            return error_json(404, "Unknown upload_id: " + upload_id);
        path     = it->second.path;
        filename = it->second.filename;
    }

    static const std::set<std::string> valid_fmts = {
        "stl", "stl_ascii", "obj", "ply", "step"
    };
    if (!valid_fmts.count(format))
        return error_json(400, "Unknown format '" + format
                        + "'. Use: stl stl_ascii obj ply step");

    try {
        GeometryEngine geo([](const std::string&, double) {});
        SurfaceManifold sm = geo.importCAD(path);
        (void)geo.repairManifold(sm);

        std::string ext      = (format == "stl_ascii") ? "stl" : format;
        std::string tmp_path = "/tmp/" + gen_id() + "_export." + ext;

        if (format == "obj") {
            auto groups = geo.identifyFlatFaces(sm);
            cfd::makeExporter("obj", &groups)->write(sm, tmp_path);
        } else {
            cfd::makeExporter(format)->write(sm, tmp_path);
        }

        std::ifstream f(tmp_path, std::ios::binary);
        if (!f.is_open())
            return error_json(500, "Export write failed");
        std::string body((std::istreambuf_iterator<char>(f)),
                          std::istreambuf_iterator<char>());
        f.close();
        ::unlink(tmp_path.c_str());

        size_t dot = filename.rfind('.');
        std::string stem = (dot == std::string::npos) ? filename : filename.substr(0, dot);
        std::string dl_name = stem + "." + ext;

        const std::string mime =
            (ext == "stl")  ? "model/stl" :
            (ext == "obj")  ? "text/plain" :
            (ext == "step") ? "application/step" :
                              "application/octet-stream";

        std::ostringstream hdr;
        hdr << "HTTP/1.1 200 OK\r\n"
            << "Content-Type: " << mime << "\r\n"
            << "Content-Disposition: attachment; filename=\""
            << json_escape(dl_name) << "\"\r\n"
            << "Content-Length: " << body.size() << "\r\n"
            << "Connection: close\r\n"
            << "Access-Control-Allow-Origin: *\r\n"
            << "\r\n"
            << body;
        return hdr.str();
    } catch (const std::exception& e) {
        return error_json(500, std::string("Export failed: ") + e.what());
    }
}

// GET /api/uploads
// Returns a JSON array of all currently stored uploads.
static std::string handle_list_uploads() {
    std::lock_guard<std::mutex> lk(store_mu);
    std::ostringstream out;
    out << "{\"uploads\":[";
    bool first = true;
    for (const auto& kv : uploads) {
        if (!first) out << ",";
        out << "{\"id\":\"" << json_escape(kv.first)
            << "\",\"filename\":\"" << json_escape(kv.second.filename)
            << "\",\"triangle_count\":" << kv.second.tri_count << "}";
        first = false;
    }
    out << "]}";
    return ok_json(out.str());
}

// DELETE /api/upload/<id>
// Removes the upload record by ID, deletes the temp file from disk.
static std::string handle_delete_upload(const std::string& upload_id) {
    std::string path;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(upload_id);
        if (it == uploads.end())
            return error_json(404, "Unknown upload_id: " + upload_id);
        path = it->second.path;
        ::unlink(path.c_str());
        uploads.erase(it);
    }
    return ok_json("{\"ok\":true}");
}

// =============================================================================
// Static file serving
// =============================================================================

static std::string WEBUI_ROOT = "webui"; // resolved at runtime relative to exe location

static std::string handle_static(const std::string& req_path) {
    // Reject path traversal attempts
    if (req_path.find("..") != std::string::npos)
        return error_json(400, "Invalid path");

    std::string rel = strip_query(req_path);
    if (!rel.empty() && rel.front() == '/') rel = rel.substr(1);

    // Root → landing page; /app → CFD app
    std::string file_path;
    if (rel.empty() || rel == "index.html") {
        file_path = WEBUI_ROOT + "/landing.html";
    } else if (rel == "app" || rel == "app/") {
        file_path = WEBUI_ROOT + "/index.html";
    } else {
        file_path = WEBUI_ROOT + "/" + rel;
    }

    struct stat st{};
    if (stat(file_path.c_str(), &st) != 0 || !S_ISREG(st.st_mode)) {
        // Fall back: unknown paths go to landing, /app/* go to CFD app
        if (rel.rfind("app", 0) == 0)
            file_path = WEBUI_ROOT + "/index.html";
        else
            file_path = WEBUI_ROOT + "/landing.html";
    }

    return serve_file(file_path, mime_for(file_path));
}

// =============================================================================
// Primitive / Boolean handlers
// =============================================================================

// Register an in-memory STL blob as an upload and return JSON {file_id, triangle_count}
static std::string register_stl_blob(const std::string& stl_data, const std::string& filename) {
    std::string id = gen_id();
    std::string tmp_path = "/tmp/" + id + ".stl";
    {
        std::ofstream out(tmp_path, std::ios::binary);
        if (!out.is_open())
            return error_json(500, "Cannot write to " + tmp_path);
        out.write(stl_data.data(), static_cast<std::streamsize>(stl_data.size()));
    }
    int tri_count = count_stl_triangles(stl_data);
    {
        std::lock_guard<std::mutex> lk(store_mu);
        uploads[id] = UploadRecord{ tmp_path, filename, tri_count, {}, {} };
    }

    // Mirror to Supabase Storage and database if configured
    if (SupabaseClient::isConfigured()) {
        std::string storage_path = id + "/" + filename;
        SupabaseClient::uploadFile("stl-files", storage_path, stl_data, "model/stl");
        std::ostringstream row;
        row << "{\"id\":\"" << id << "\",\"filename\":\"" << json_escape(filename)
            << "\",\"tri_count\":" << tri_count << "}";
        SupabaseClient::insertRow("uploads", row.str());
    }

    std::ostringstream resp;
    resp << "{\"file_id\":\"" << id
         << "\",\"triangle_count\":" << tri_count
         << ",\"filename\":\"" << json_escape(filename) << "\"}";
    return ok_json(resp.str());
}

// POST /api/primitives/box      body: {cx,cy,cz, hx,hy,hz}
// POST /api/primitives/sphere   body: {cx,cy,cz, rx,ry,rz, res}
// POST /api/primitives/cylinder body: {cx,cy,cz, radius, height, segs}
// POST /api/primitives/cone     body: {cx,cy,cz, radius, height, segs}
// POST /api/primitives/torus    body: {cx,cy,cz, major_r, minor_r, major_segs, minor_segs}
// POST /api/primitives/plane    body: {cx,cy,cz, width, depth, subdiv}
static std::string handle_primitive(const std::string& type, const HttpRequest& req) {
    if (type == "box") {
        double cx = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double hx = json_double_value(req.body,"hx",0.5), hy = json_double_value(req.body,"hy",0.5), hz = json_double_value(req.body,"hz",0.5);
        return register_stl_blob(generate_box_stl(cx,cy,cz,hx,hy,hz), "box.stl");
    }
    if (type == "sphere") {
        double cx = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double rx = json_double_value(req.body,"rx",0.5), ry = json_double_value(req.body,"ry",0.5), rz = json_double_value(req.body,"rz",0.5);
        int    res = json_int_value(req.body,"res",32);
        return register_stl_blob(generate_ellipsoid_stl(cx,cy,cz,rx,ry,rz,res), "sphere.stl");
    }
    if (type == "cylinder") {
        double cx     = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double radius = json_double_value(req.body,"radius",0.5);
        double height = json_double_value(req.body,"height",1.0);
        int    segs   = json_int_value(req.body,"segs",32);
        return register_stl_blob(cylinder_stl(cx,cy,cz,radius,height,segs), "cylinder.stl");
    }
    if (type == "cone") {
        double cx     = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double radius = json_double_value(req.body,"radius",0.5);
        double height = json_double_value(req.body,"height",1.0);
        int    segs   = json_int_value(req.body,"segs",24);
        return register_stl_blob(cone_stl(cx,cy,cz,radius,height,segs), "cone.stl");
    }
    if (type == "torus") {
        double cx      = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double major_r = json_double_value(req.body,"major_r",1.0);
        double minor_r = json_double_value(req.body,"minor_r",0.25);
        int major_segs = json_int_value(req.body,"major_segs",32);
        int minor_segs = json_int_value(req.body,"minor_segs",16);
        return register_stl_blob(torus_stl(cx,cy,cz,major_r,minor_r,major_segs,minor_segs), "torus.stl");
    }
    if (type == "plane") {
        double cx    = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double width = json_double_value(req.body,"width",1.0);
        double depth = json_double_value(req.body,"depth",1.0);
        int    subdiv = json_int_value(req.body,"subdiv",4);
        return register_stl_blob(plane_stl(cx,cy,cz,width,depth,subdiv), "plane.stl");
    }
    if (type == "icosphere") {
        double cx     = json_double_value(req.body,"cx",0), cy = json_double_value(req.body,"cy",0), cz = json_double_value(req.body,"cz",0);
        double radius = json_double_value(req.body,"radius",0.5);
        int    subdiv = json_int_value(req.body,"subdivisions",2);
        return register_stl_blob(icosphere_stl(cx,cy,cz,radius,subdiv), "icosphere.stl");
    }
    return error_json(400, "Unknown primitive: " + type);
}

// =============================================================================
// v1 API — arbitrary-mesh boolean jobs, bifilar jobs, pipeline triggers
// (Lovable + Supabase compatible: async job pattern, CORS, JSON)
// =============================================================================

struct BooleanMeshJob {
    std::string status;        // "pending"|"running"|"complete"|"error"
    std::string message;
    double      progress = 0.0;
    std::string mesh_a_id;
    std::string mesh_b_id;
    std::string result_stl_path;
    std::string result_vtu_path;
    std::string supabase_url;
    std::mutex  mu;
};

struct BifilarJobV1 {
    std::string status;
    std::string message;
    double      progress = 0.0;
    double R=1.0, H=4.0, revolutions=2.0;
    int    resolution=40;
    std::string result_stl_path;
    std::string result_vtu_path;
    std::string supabase_url;
    std::mutex  mu;
};

struct BifilarVisualJob {
    double R=1.0, H=4.0, revolutions=2.0;
    int    resolution=40;
    std::string status, message;
    double      progress = 0.0;
    std::string inner1_path, inner2_path;
    std::mutex  mu;
};

static std::map<std::string, std::shared_ptr<BooleanMeshJob>>  v1_bool_jobs;
static std::map<std::string, std::shared_ptr<BifilarJobV1>>    v1_bif_jobs;
static std::map<std::string, std::shared_ptr<BifilarVisualJob>> v1_bif_vis_jobs;
static std::mutex v1_bool_mu, v1_bif_mu, v1_bif_vis_mu;

static void ensure_dir_sv(const std::string& d){
    std::string cmd = "mkdir -p \"" + d + "\"";
    std::system(cmd.c_str());
}

static std::string job_status_json(const std::string& id, const std::string& status,
                                    double progress, const std::string& message,
                                    const std::string& stl_url = "",
                                    const std::string& vtu_url = "",
                                    const std::string& sb_url  = "") {
    std::ostringstream o;
    o << "{\"job_id\":\"" << json_escape(id) << "\""
      << ",\"status\":\"" << json_escape(status) << "\""
      << ",\"progress\":" << progress
      << ",\"message\":\"" << json_escape(message) << "\"";
    if(!stl_url.empty())
        o << ",\"result\":{\"stl_url\":\"" << json_escape(stl_url) << "\""
          << ",\"vtu_url\":\"" << json_escape(vtu_url) << "\""
          << ",\"supabase_url\":\"" << json_escape(sb_url) << "\"}";
    o << "}";
    return o.str();
}

// POST /api/v1/boolean/subtract-meshes
// body: {"mesh_a_id":"...", "mesh_b_id":"..."}
static std::string handle_v1_boolean_subtract(const HttpRequest& req) {
    std::string aid = json_string_value(req.body, "mesh_a_id");
    std::string bid = json_string_value(req.body, "mesh_b_id");
    if(aid.empty() || bid.empty())
        return error_json(400, "mesh_a_id and mesh_b_id required");

    std::string path_a, path_b;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto ia = uploads.find(aid), ib = uploads.find(bid);
        if(ia == uploads.end()) return error_json(404, "mesh_a_id not found");
        if(ib == uploads.end()) return error_json(404, "mesh_b_id not found");
        path_a = ia->second.path;
        path_b = ib->second.path;
    }

    std::string jid = gen_id();
    auto job = std::make_shared<BooleanMeshJob>();
    job->status = "pending"; job->mesh_a_id = aid; job->mesh_b_id = bid;
    {
        std::lock_guard<std::mutex> lk(v1_bool_mu);
        v1_bool_jobs[jid] = job;
    }

    std::thread([jid, job, path_a, path_b](){
        auto set = [&](const std::string& s, double p, const std::string& m=""){
            std::lock_guard<std::mutex> lk(job->mu);
            job->status=s; job->progress=p; job->message=m;
        };
        try {
            set("running", 0.05, "loading meshes");
            GeometryEngine geo;
            auto ma = geo.importCAD(path_a); geo.repairManifold(ma);
            auto mb = geo.importCAD(path_b); geo.repairManifold(mb);

            set("running", 0.25, "building SDF");
            MeshSDF sdf_a(ma), sdf_b(mb);
            Vec3d la,ha,lb,hb;
            sdf_a.bounds(la,ha); sdf_b.bounds(lb,hb);
            double lo[3]={std::min(la[0],lb[0]),std::min(la[1],lb[1]),std::min(la[2],lb[2])};
            double hi[3]={std::max(ha[0],hb[0]),std::max(ha[1],hb[1]),std::max(ha[2],hb[2])};

            set("running", 0.40, "marching cubes");
            auto sdf=[&](double x,double y,double z){ return std::max(sdf_a.eval({x,y,z}),-sdf_b.eval({x,y,z})); };
            auto tris = marching_cubes_run(sdf,lo[0],lo[1],lo[2],hi[0],hi[1],hi[2],40);
            auto result = tris_to_manifold(tris);
            geo.mergeVertices(result,1e-8); geo.repairManifold(result);

            set("running", 0.65, "tetrahedral mesh");
            MeshEngine eng;
            BLParameters bl; bl.first_cell_height=0.02; bl.growth_ratio=1.2; bl.num_layers=5;
            bl.total_thickness=bl.first_cell_height*(std::pow(bl.growth_ratio,bl.num_layers)-1)/(bl.growth_ratio-1);
            auto prisms = eng.inflatePrisms(result, bl);
            auto vol    = eng.generateTetCore(prisms, 0.08);

            set("running", 0.85, "exporting");
            std::string dir = "output/boolean/" + jid;
            ensure_dir_sv(dir);
            std::string stl_p = dir+"/result.stl", vtu_p = dir+"/result.vtu";
            makeExporter("stl_binary")->write(result, stl_p);
            VTUExporter vtu_exp; vtu_exp.writeVolume(vol, vtu_p);

            std::string sb_url;
            if(SupabaseClient::isConfigured()){
                std::ifstream vf(vtu_p, std::ios::binary);
                std::string data((std::istreambuf_iterator<char>(vf)),{});
                sb_url = SupabaseClient::uploadFile("meshes", jid+"/result.vtu", data, "application/octet-stream");
            }

            std::lock_guard<std::mutex> lk(job->mu);
            job->status="complete"; job->progress=1.0;
            job->result_stl_path=stl_p; job->result_vtu_path=vtu_p;
            job->supabase_url=sb_url;
        } catch(const std::exception& e){
            set("error",0.0,e.what());
        }
    }).detach();

    return accepted_json("{\"job_id\":\""+jid+"\",\"status\":\"pending\"}");
}

// GET /api/v1/boolean/jobs/:id
static std::string handle_v1_boolean_status(const std::string& jid) {
    std::shared_ptr<BooleanMeshJob> job;
    { std::lock_guard<std::mutex> lk(v1_bool_mu); auto it=v1_bool_jobs.find(jid); if(it==v1_bool_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    std::string stl_url = job->status=="complete" ? "/api/v1/boolean/jobs/"+jid+"/result.stl" : "";
    std::string vtu_url = job->status=="complete" ? "/api/v1/boolean/jobs/"+jid+"/result.vtu" : "";
    return ok_json(job_status_json(jid,job->status,job->progress,job->message,stl_url,vtu_url,job->supabase_url));
}

// GET /api/v1/boolean/jobs/:id/result.stl
static std::string handle_v1_boolean_stl(const std::string& jid) {
    std::shared_ptr<BooleanMeshJob> job;
    { std::lock_guard<std::mutex> lk(v1_bool_mu); auto it=v1_bool_jobs.find(jid); if(it==v1_bool_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    if(job->status!="complete") return error_json(409,"job not complete");
    return serve_file(job->result_stl_path, "model/stl");
}

// GET /api/v1/boolean/jobs/:id/result.vtu
static std::string handle_v1_boolean_vtu(const std::string& jid) {
    std::shared_ptr<BooleanMeshJob> job;
    { std::lock_guard<std::mutex> lk(v1_bool_mu); auto it=v1_bool_jobs.find(jid); if(it==v1_bool_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    if(job->status!="complete") return error_json(409,"job not complete");
    return serve_file(job->result_vtu_path, "application/xml");
}

// POST /api/v1/bifilar/generate
// body: {"R":1.0, "H":4.0, "revolutions":2.0, "resolution":40}
static std::string handle_v1_bifilar_generate(const HttpRequest& req) {
    auto job = std::make_shared<BifilarJobV1>();
    job->R           = json_double_value(req.body,"R",1.0);
    job->H           = json_double_value(req.body,"H",4.0);
    job->revolutions = json_double_value(req.body,"revolutions",2.0);
    job->resolution  = json_int_value(req.body,"resolution",40);
    job->status = "pending";

    std::string jid = gen_id();
    { std::lock_guard<std::mutex> lk(v1_bif_mu); v1_bif_jobs[jid] = job; }

    std::thread([jid, job](){
        auto set=[&](const std::string& s, double p, const std::string& m=""){
            std::lock_guard<std::mutex> lk(job->mu);
            job->status=s; job->progress=p; job->message=m;
        };
        try {
            set("running",0.05,"generating geometry");
            BifilarConfig cfg;
            cfg.R=job->R; cfg.H=job->H; cfg.revolutions=job->revolutions; cfg.resolution=job->resolution;
            auto surface = generateBifilarSurface(cfg);

            set("running",0.35,"boundary layer inflation");
            MeshEngine eng;
            BLParameters bl;
            bl.first_cell_height=0.008*cfg.R; bl.growth_ratio=1.1; bl.num_layers=10;
            bl.total_thickness=bl.first_cell_height*(std::pow(bl.growth_ratio,10)-1)/(bl.growth_ratio-1);
            auto prisms = eng.inflatePrisms(surface,bl);

            set("running",0.65,"tet core");
            auto vol = eng.generateTetCore(prisms, 0.06*cfg.R);

            set("running",0.85,"exporting");
            std::string dir="output/bifilar/"+jid;
            ensure_dir_sv(dir);
            std::string stl_p=dir+"/surface.stl", vtu_p=dir+"/mesh.vtu";
            makeExporter("stl_binary")->write(surface, stl_p);
            VTUExporter vtu; vtu.writeVolume(vol, vtu_p);

            std::string sb_url;
            if(SupabaseClient::isConfigured()){
                std::ifstream vf(vtu_p,std::ios::binary);
                std::string data((std::istreambuf_iterator<char>(vf)),{});
                sb_url=SupabaseClient::uploadFile("meshes",jid+"/bifilar.vtu",data,"application/octet-stream");
            }

            std::lock_guard<std::mutex> lk(job->mu);
            job->status="complete"; job->progress=1.0;
            job->result_stl_path=stl_p; job->result_vtu_path=vtu_p;
            job->supabase_url=sb_url;
        } catch(const std::exception& e){
            set("error",0.0,e.what());
        }
    }).detach();

    return accepted_json("{\"job_id\":\""+jid+"\",\"status\":\"pending\"}");
}

// GET /api/v1/bifilar/jobs/:id
static std::string handle_v1_bifilar_status(const std::string& jid) {
    std::shared_ptr<BifilarJobV1> job;
    { std::lock_guard<std::mutex> lk(v1_bif_mu); auto it=v1_bif_jobs.find(jid); if(it==v1_bif_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    std::string stl_url = job->status=="complete" ? "/api/v1/bifilar/jobs/"+jid+"/result.stl" : "";
    std::string vtu_url = job->status=="complete" ? "/api/v1/bifilar/jobs/"+jid+"/result.vtu" : "";
    return ok_json(job_status_json(jid,job->status,job->progress,job->message,stl_url,vtu_url,job->supabase_url));
}

// GET /api/v1/bifilar/jobs/:id/result.stl
static std::string handle_v1_bifilar_stl(const std::string& jid) {
    std::shared_ptr<BifilarJobV1> job;
    { std::lock_guard<std::mutex> lk(v1_bif_mu); auto it=v1_bif_jobs.find(jid); if(it==v1_bif_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    if(job->status!="complete") return error_json(409,"job not complete");
    return serve_file(job->result_stl_path,"model/stl");
}

// GET /api/v1/bifilar/jobs/:id/result.vtu
static std::string handle_v1_bifilar_vtu(const std::string& jid) {
    std::shared_ptr<BifilarJobV1> job;
    { std::lock_guard<std::mutex> lk(v1_bif_mu); auto it=v1_bif_jobs.find(jid); if(it==v1_bif_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    if(job->status!="complete") return error_json(409,"job not complete");
    return serve_file(job->result_vtu_path,"application/xml");
}

// POST /api/v1/bifilar/generate-visual
// Generates the two solid inner helical cylinders for colour visualisation.
static std::string handle_v1_bifilar_visual(const HttpRequest& req) {
    auto job = std::make_shared<BifilarVisualJob>();
    job->R           = json_double_value(req.body,"R",1.0);
    job->H           = json_double_value(req.body,"H",4.0);
    job->revolutions = json_double_value(req.body,"revolutions",2.0);
    job->resolution  = json_int_value(req.body,"resolution",40);
    job->status = "pending";

    std::string jid = gen_id();
    { std::lock_guard<std::mutex> lk(v1_bif_vis_mu); v1_bif_vis_jobs[jid] = job; }

    std::thread([jid, job](){
        auto set = [&](const std::string& s, double p, const std::string& m=""){
            std::lock_guard<std::mutex> lk(job->mu);
            job->status=s; job->progress=p; job->message=m;
        };
        try {
            BifilarConfig cfg;
            cfg.R=job->R; cfg.H=job->H; cfg.revolutions=job->revolutions; cfg.resolution=job->resolution;

            set("running",0.05,"generating inner cylinder 1");
            auto surf1 = generateBifilarInnerCylinder(cfg,1);

            set("running",0.55,"generating inner cylinder 2");
            auto surf2 = generateBifilarInnerCylinder(cfg,2);

            set("running",0.90,"exporting STLs");
            std::string dir="output/bifilar/"+jid;
            ensure_dir_sv(dir);
            std::string p1=dir+"/inner1.stl", p2=dir+"/inner2.stl";
            makeExporter("stl_binary")->write(surf1,p1);
            makeExporter("stl_binary")->write(surf2,p2);

            std::lock_guard<std::mutex> lk(job->mu);
            job->status="complete"; job->progress=1.0;
            job->inner1_path=p1; job->inner2_path=p2;
        } catch(const std::exception& e){
            set("error",0.0,e.what());
        }
    }).detach();

    return accepted_json("{\"job_id\":\""+jid+"\",\"status\":\"pending\"}");
}

// GET /api/v1/bifilar/visual-jobs/:id
static std::string handle_v1_bifilar_visual_status(const std::string& jid) {
    std::shared_ptr<BifilarVisualJob> job;
    { std::lock_guard<std::mutex> lk(v1_bif_vis_mu); auto it=v1_bif_vis_jobs.find(jid); if(it==v1_bif_vis_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    std::ostringstream o;
    o << "{\"job_id\":\"" << jid << "\",\"status\":\"" << job->status << "\""
      << ",\"progress\":" << job->progress << ",\"message\":\"" << json_escape(job->message) << "\"";
    if (job->status=="complete")
        o << ",\"inner1_url\":\"/api/v1/bifilar/visual-jobs/"+jid+"/inner1.stl\""
          << ",\"inner2_url\":\"/api/v1/bifilar/visual-jobs/"+jid+"/inner2.stl\"";
    o << "}";
    return ok_json(o.str());
}

// GET /api/v1/bifilar/visual-jobs/:id/inner1.stl  or  inner2.stl
static std::string handle_v1_bifilar_visual_stl(const std::string& jid, const std::string& which) {
    std::shared_ptr<BifilarVisualJob> job;
    { std::lock_guard<std::mutex> lk(v1_bif_vis_mu); auto it=v1_bif_vis_jobs.find(jid); if(it==v1_bif_vis_jobs.end()) return error_json(404,"job not found"); job=it->second; }
    std::lock_guard<std::mutex> lk(job->mu);
    if(job->status!="complete") return error_json(409,"job not complete");
    const std::string& path = (which=="inner1.stl") ? job->inner1_path : job->inner2_path;
    return serve_file(path,"model/stl");
}

// GET /api/v1/task1/models
static std::string handle_v1_task1_models() {
    const std::string BASE="https://huggingface.co/datasets/Thingi10K/Thingi10K/resolve/main/raw_meshes/";
    struct M { const char* id; const char* name; };
    M ms[] = {{"100035","Tail Section"},{"100075","Robot Arm Front"},
               {"100335","Stonehenge Stone"},{"100388","Menou Final"},
               {"1005587","Bird Feeder"}};
    std::ostringstream o; o << "{\"models\":[";
    for(int i=0;i<5;i++){
        if(i) o<<",";
        std::string path="output/task1/models/"+std::string(ms[i].id)+".stl";
        struct stat st{}; bool av=(stat(path.c_str(),&st)==0);
        o << "{\"id\":\"" << ms[i].id << "\",\"name\":\"" << ms[i].name << "\""
          << ",\"path\":\"" << json_escape(path) << "\",\"available\":" << (av?"true":"false") << "}";
    }
    o << "]}";
    return ok_json(o.str());
}

// POST /api/v1/task1/run
static std::string handle_v1_task1_run(const HttpRequest&) {
    std::string jid = gen_id();
    std::thread([jid](){
        std::string cmd = "./build/task1_pipeline --output-dir output/task1 > output/task1/run_"+jid+".log 2>&1";
        ensure_dir_sv("output/task1");
        std::system(cmd.c_str());
    }).detach();
    return accepted_json("{\"job_id\":\""+jid+"\",\"status\":\"running\","
                         "\"log\":\"output/task1/run_"+jid+".log\"}");
}

// POST /api/v1/task2/run
static std::string handle_v1_task2_run(const HttpRequest& req) {
    double R=json_double_value(req.body,"R",1.0);
    double H=json_double_value(req.body,"H",4.0);
    double revs=json_double_value(req.body,"revolutions",2.0);
    int res=json_int_value(req.body,"resolution",40);
    std::string jid = gen_id();
    std::thread([jid,R,H,revs,res](){
        ensure_dir_sv("output/task2");
        std::string cmd = "./build/task2_pipeline --output-dir output/task2"
            " --R "+std::to_string(R)+" --H "+std::to_string(H)
            +" --revs "+std::to_string(revs)+" --resolution "+std::to_string(res)
            +" > output/task2/run_"+jid+".log 2>&1";
        std::system(cmd.c_str());
    }).detach();
    return accepted_json("{\"job_id\":\""+jid+"\",\"status\":\"running\","
                         "\"log\":\"output/task2/run_"+jid+".log\"}");
}

// ─────────────────────────────────────────────────────────────────────────────

// POST /api/boolean/subtract
// body: {bcx,bcy,bcz,bhx,bhy,bhz, scx,scy,scz,srx,sry,srz, resolution}
static std::string handle_boolean_subtract(const HttpRequest& req) {
    double bcx=json_double_value(req.body,"bcx",0),  bcy=json_double_value(req.body,"bcy",0),  bcz=json_double_value(req.body,"bcz",0);
    double bhx=json_double_value(req.body,"bhx",0.5),bhy=json_double_value(req.body,"bhy",0.5),bhz=json_double_value(req.body,"bhz",0.5);
    double scx=json_double_value(req.body,"scx",0),  scy=json_double_value(req.body,"scy",0),  scz=json_double_value(req.body,"scz",0);
    double srx=json_double_value(req.body,"srx",0.5),sry=json_double_value(req.body,"sry",0.5),srz=json_double_value(req.body,"srz",0.5);
    int res = json_int_value(req.body,"resolution",60);
    std::string stl = boolean_subtract(bcx,bcy,bcz,bhx,bhy,bhz, scx,scy,scz,srx,sry,srz, res);
    return register_stl_blob(stl, "boolean_subtract.stl");
}

// =============================================================================
// Render handlers
// =============================================================================

// Parse a hex color string like "#1a1a2e" or "1a1a2e" into 0xRRGGBB
static uint32_t parse_hex_color(const std::string& s, uint32_t def) {
    std::string h = s;
    if (!h.empty() && h[0] == '#') h = h.substr(1);
    if (h.size() != 6) return def;
    try {
        return static_cast<uint32_t>(std::stoul(h, nullptr, 16));
    } catch (...) { return def; }
}

// Parse a JSON array like [1.0, 2.0, 3.0] at a given key into out[3]
// Returns true if exactly 3 values found.
static bool json_vec3(const std::string& json, const std::string& key, double out[3]) {
    std::string search = "\"" + key + "\"";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return false;
    pos += search.size();
    while (pos < json.size() && (json[pos]==' '||json[pos]==':'||json[pos]=='\t')) ++pos;
    if (pos >= json.size() || json[pos] != '[') return false;
    ++pos;
    for (int i = 0; i < 3; ++i) {
        while (pos < json.size() && (json[pos]==' '||json[pos]==','||json[pos]=='\t'||json[pos]=='\n')) ++pos;
        if (pos >= json.size()) return false;
        try {
            size_t consumed = 0;
            out[i] = std::stod(json.substr(pos), &consumed);
            pos += consumed;
        } catch (...) { return false; }
    }
    return true;
}

// Read binary STL from file path and return flat triangle array (N*9 floats)
static bool load_stl_triangles(const std::string& path, std::vector<float>& tris) {
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return false;
    f.seekg(0, std::ios::end);
    size_t sz = static_cast<size_t>(f.tellg());
    f.seekg(0, std::ios::beg);
    if (sz < 84) return false;
    std::vector<char> buf(sz);
    f.read(buf.data(), static_cast<std::streamsize>(sz));
    if (!f) return false;

    // Binary STL: 80-byte header, 4-byte count, then 50 bytes per tri
    uint32_t count = 0;
    std::memcpy(&count, buf.data() + 80, 4);
    if (84 + static_cast<size_t>(count) * 50 > sz) {
        // Could be ASCII; try to parse facets
        std::string text(buf.begin(), buf.end());
        tris.clear();
        size_t p = 0;
        auto skip_to = [&](const std::string& kw) -> bool {
            size_t q = text.find(kw, p);
            if (q == std::string::npos) return false;
            p = q + kw.size();
            return true;
        };
        while (skip_to("vertex")) {
            // read 3 floats
            for (int v = 0; v < 3; ++v) {
                if (v > 0 && !skip_to("vertex")) break;
                std::istringstream ss(text.substr(p, 80));
                float x, y, z;
                if (!(ss >> x >> y >> z)) break;
                tris.push_back(x); tris.push_back(y); tris.push_back(z);
                p += static_cast<size_t>(ss.tellg()) <= 0 ? 0 : static_cast<size_t>(ss.tellg());
            }
        }
        return !tris.empty();
    }

    tris.clear();
    tris.reserve(count * 9);
    for (uint32_t i = 0; i < count; ++i) {
        const char* tri = buf.data() + 84 + i * 50;
        // skip 12-byte normal, then 3 verts × 12 bytes
        for (int v = 0; v < 3; ++v) {
            float x, y, z;
            std::memcpy(&x, tri + 12 + v*12 + 0, 4);
            std::memcpy(&y, tri + 12 + v*12 + 4, 4);
            std::memcpy(&z, tri + 12 + v*12 + 8, 4);
            tris.push_back(x); tris.push_back(y); tris.push_back(z);
        }
    }
    return !tris.empty();
}

// Base64 decode helper (for stl_data field)
static std::string base64_decode(const std::string& in) {
    static const std::string chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    auto is_b64 = [](unsigned char c) {
        return (isalnum(c) || c == '+' || c == '/');
    };
    std::string out;
    int i = 0;
    unsigned char buf4[4], buf3[3];
    size_t pos = 0;
    while (pos < in.size() && in[pos] != '=' && is_b64((unsigned char)in[pos])) {
        buf4[i++] = (unsigned char)in[pos++];
        if (i == 4) {
            for (int j = 0; j < 4; ++j)
                buf4[j] = static_cast<unsigned char>(chars.find((char)buf4[j]));
            buf3[0] = (buf4[0] << 2) | (buf4[1] >> 4);
            buf3[1] = ((buf4[1] & 0xF) << 4) | (buf4[2] >> 2);
            buf3[2] = ((buf4[2] & 0x3) << 6) | buf4[3];
            for (int j = 0; j < 3; ++j) out += (char)buf3[j];
            i = 0;
        }
    }
    if (i > 0) {
        for (int j = i; j < 4; ++j) buf4[j] = 0;
        for (int j = 0; j < 4; ++j)
            buf4[j] = static_cast<unsigned char>(chars.find((char)buf4[j]) == std::string::npos ? 0 : chars.find((char)buf4[j]));
        buf3[0] = (buf4[0] << 2) | (buf4[1] >> 4);
        buf3[1] = ((buf4[1] & 0xF) << 4) | (buf4[2] >> 2);
        buf3[2] = ((buf4[2] & 0x3) << 6) | buf4[3];
        for (int j = 0; j < i - 1; ++j) out += (char)buf3[j];
    }
    return out;
}

// GET /api/render/defaults
static std::string handle_render_defaults() {
    return ok_json(R"({"eye":[2,2,2],"target":[0,0,0],"up":[0,1,0],"fov":45,"width":800,"height":600})");
}

// POST /api/render
static std::string handle_render(const HttpRequest& req) {
    // Parse camera
    render::Camera cam;
    cam.eye[0]=2; cam.eye[1]=2; cam.eye[2]=2;
    cam.target[0]=cam.target[1]=cam.target[2]=0;
    cam.up[0]=0; cam.up[1]=1; cam.up[2]=0;
    cam.fov_deg = 45.0;
    cam.width = 800;
    cam.height = 600;

    double tmp3[3];
    if (json_vec3(req.body, "eye",    tmp3)) { cam.eye[0]=tmp3[0]; cam.eye[1]=tmp3[1]; cam.eye[2]=tmp3[2]; }
    if (json_vec3(req.body, "target", tmp3)) { cam.target[0]=tmp3[0]; cam.target[1]=tmp3[1]; cam.target[2]=tmp3[2]; }
    if (json_vec3(req.body, "up",     tmp3)) { cam.up[0]=tmp3[0]; cam.up[1]=tmp3[1]; cam.up[2]=tmp3[2]; }

    double fov = json_double_value(req.body, "fov", 45.0);
    if (fov > 0) cam.fov_deg = fov;
    int w = json_int_value(req.body, "width",  800);
    int h = json_int_value(req.body, "height", 600);
    if (w > 0 && w <= 4096) cam.width  = w;
    if (h > 0 && h <= 4096) cam.height = h;

    uint32_t bg_col   = 0x1a1a2e;
    uint32_t mesh_col = 0x6b58cc;
    {
        std::string bgs = json_string_value(req.body, "bg_color");
        if (!bgs.empty()) bg_col = parse_hex_color(bgs, bg_col);
        std::string mcs = json_string_value(req.body, "mesh_color");
        if (!mcs.empty()) mesh_col = parse_hex_color(mcs, mesh_col);
    }

    // Load geometry
    std::vector<float> tris;

    std::string upload_id = json_string_value(req.body, "upload_id");
    if (!upload_id.empty()) {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(upload_id);
        if (it == uploads.end())
            return error_json(404, "Unknown upload_id: " + upload_id);
        std::string path = it->second.path;
        lk.~lock_guard();
        // Load STL from file
        if (!load_stl_triangles(path, tris))
            return error_json(500, "Failed to load STL from upload");
    } else {
        // Try base64 stl_data
        std::string stl_b64 = json_string_value(req.body, "stl_data");
        if (!stl_b64.empty()) {
            std::string stl_bytes = base64_decode(stl_b64);
            // Write to temp file and load
            std::string tmp_path = "/tmp/render_" + gen_id() + ".stl";
            {
                std::ofstream tf(tmp_path, std::ios::binary);
                tf.write(stl_bytes.data(), static_cast<std::streamsize>(stl_bytes.size()));
            }
            if (!load_stl_triangles(tmp_path, tris))
                return error_json(400, "Failed to decode stl_data");
            ::unlink(tmp_path.c_str());
        }
    }

    // Render
    auto png = render::render_stl(tris, cam, bg_col, mesh_col);

    // Return raw PNG with appropriate headers
    std::string body(reinterpret_cast<const char*>(png.data()), png.size());
    std::ostringstream resp;
    resp << "HTTP/1.1 200 OK\r\n"
         << "Content-Type: image/png\r\n"
         << "Content-Length: " << png.size() << "\r\n"
         << "Connection: close\r\n"
         << "Access-Control-Allow-Origin: *\r\n"
         << "\r\n";
    std::string header = resp.str();
    return header + body;
}

// =============================================================================
// Node shader compiler handler
// =============================================================================

// POST /api/shaders/from_nodes
// Body: node graph JSON (see node_compiler.hpp for format)
// Returns: {"vertex":"...","fragment":"...","uniforms":{...},"error":""}
static std::string handle_compile_nodes(const HttpRequest& req) {
    if (req.body.empty())
        return error_json(400, "Request body is empty");

    std::string result = node_compiler::compile(req.body);

    // Persist to Supabase if configured
    if (SupabaseClient::isConfigured()) {
        // Build a flat row; store the full compiled JSON as compiled_json
        std::ostringstream row;
        row << "{\"graph_json\":\""    << json_escape(req.body)  << "\""
            << ",\"compiled_json\":\"" << json_escape(result)    << "\"}";
        SupabaseClient::insertRow("node_shaders", row.str());
    }

    return http_response(200, "OK", "application/json", result);
}

// =============================================================================
// Render frame handler (POST /api/render/frame)
// =============================================================================

// Apply a 3-D Euler rotation (Rx * Ry * Rz) + translation + scale to a point.
static void transform_point(const double p_in[3], double p_out[3],
                             const double rx, const double ry, const double rz,
                             const double tx, const double ty, const double tz,
                             const double sx, const double sy, const double sz) {
    // Rotation around X
    double cx = std::cos(rx), sx_ = std::sin(rx);
    double y1 = cx * p_in[1] - sx_ * p_in[2];
    double z1 = sx_ * p_in[1] + cx * p_in[2];

    // Rotation around Y
    double cy = std::cos(ry), sy_ = std::sin(ry);
    double x2 = cy * p_in[0] + sy_ * z1;
    double z2 = -sy_ * p_in[0] + cy * z1;

    // Rotation around Z
    double cz = std::cos(rz), sz_ = std::sin(rz);
    double x3 = cz * x2 - sz_ * y1;
    double y3 = sz_ * x2 + cz * y1;
    double z3 = z2;

    // Translate then scale
    p_out[0] = (x3 + tx) * sx;
    p_out[1] = (y3 + ty) * sy;
    p_out[2] = (z3 + tz) * sz;
}

// Persistent frame counter per file_id
static std::map<std::string, int>  frame_counters;
static std::mutex                  frame_counter_mu;

// POST /api/render/frame
// Body JSON: {
//   "file_id":"...", "position":[x,y,z], "rotation_euler":[rx,ry,rz],
//   "scale":[sx,sy,sz], "eye":[x,y,z], "target":[x,y,z],
//   "width":800, "height":600, "bg_color":"#1a1a2e", "mesh_color":"#6b58cc"
// }
// Returns JSON: {"url":"...", "frame_number":N}  (or raw PNG if no Supabase)
static std::string handle_render_frame(const HttpRequest& req) {
    std::string file_id = json_string_value(req.body, "file_id");
    if (file_id.empty())
        return error_json(400, "file_id required");

    // Load triangles
    std::vector<float> tris;
    std::string stl_path;
    {
        std::lock_guard<std::mutex> lk(store_mu);
        auto it = uploads.find(file_id);
        if (it == uploads.end())
            return error_json(404, "Unknown file_id: " + file_id);
        stl_path = it->second.path;
    }
    if (!load_stl_triangles(stl_path, tris))
        return error_json(500, "Failed to load STL for file_id: " + file_id);

    // Parse transform params
    double pos[3]  = {0,0,0};
    double rot[3]  = {0,0,0};
    double scl[3]  = {1,1,1};
    json_vec3(req.body, "position",       pos);
    json_vec3(req.body, "rotation_euler", rot);
    json_vec3(req.body, "scale",          scl);

    // Apply transform to every triangle vertex
    std::vector<float> transformed;
    transformed.reserve(tris.size());
    for (size_t i = 0; i + 8 < tris.size() + 1; i += 9) {
        for (int v = 0; v < 3; ++v) {
            double pin[3] = { tris[i + v*3 + 0],
                              tris[i + v*3 + 1],
                              tris[i + v*3 + 2] };
            double pout[3];
            transform_point(pin, pout,
                            rot[0], rot[1], rot[2],
                            pos[0], pos[1], pos[2],
                            scl[0], scl[1], scl[2]);
            transformed.push_back(static_cast<float>(pout[0]));
            transformed.push_back(static_cast<float>(pout[1]));
            transformed.push_back(static_cast<float>(pout[2]));
        }
    }

    // Camera
    render::Camera cam;
    cam.eye[0]=2; cam.eye[1]=2; cam.eye[2]=2;
    cam.target[0]=cam.target[1]=cam.target[2]=0;
    cam.up[0]=0; cam.up[1]=1; cam.up[2]=0;
    cam.fov_deg = 45.0;
    cam.width = 800; cam.height = 600;

    double tmp3[3];
    if (json_vec3(req.body, "eye",    tmp3)) { cam.eye[0]=tmp3[0]; cam.eye[1]=tmp3[1]; cam.eye[2]=tmp3[2]; }
    if (json_vec3(req.body, "target", tmp3)) { cam.target[0]=tmp3[0]; cam.target[1]=tmp3[1]; cam.target[2]=tmp3[2]; }

    int w = json_int_value(req.body, "width",  800);
    int h = json_int_value(req.body, "height", 600);
    if (w > 0 && w <= 4096) cam.width  = w;
    if (h > 0 && h <= 4096) cam.height = h;

    uint32_t bg_col   = 0x1a1a2e;
    uint32_t mesh_col = 0x6b58cc;
    {
        std::string bgs = json_string_value(req.body, "bg_color");
        if (!bgs.empty()) bg_col = parse_hex_color(bgs, bg_col);
        std::string mcs = json_string_value(req.body, "mesh_color");
        if (!mcs.empty()) mesh_col = parse_hex_color(mcs, mesh_col);
    }

    // Render
    auto png = render::render_stl(transformed, cam, bg_col, mesh_col);

    // Frame counter
    int frame_n = 0;
    {
        std::lock_guard<std::mutex> lk(frame_counter_mu);
        frame_n = ++frame_counters[file_id];
    }

    // If Supabase configured: upload PNG and return JSON with URL
    if (SupabaseClient::isConfigured()) {
        std::string storage_path = file_id + "/frame_" + std::to_string(frame_n) + ".png";
        std::string png_data(reinterpret_cast<const char*>(png.data()), png.size());
        std::string url = SupabaseClient::uploadFile("renders", storage_path, png_data, "image/png");
        if (url.empty())
            url = SupabaseClient::getPublicUrl("renders", storage_path);

        std::ostringstream resp_body;
        resp_body << "{\"url\":\"" << json_escape(url) << "\""
                  << ",\"frame_number\":" << frame_n << "}";
        return http_response(200, "OK", "application/json", resp_body.str());
    }

    // Otherwise return raw PNG
    std::string body(reinterpret_cast<const char*>(png.data()), png.size());
    std::ostringstream resp;
    resp << "HTTP/1.1 200 OK\r\n"
         << "Content-Type: image/png\r\n"
         << "Content-Length: " << png.size() << "\r\n"
         << "Connection: close\r\n"
         << "Access-Control-Allow-Origin: *\r\n"
         << "X-Frame-Number: " << frame_n << "\r\n"
         << "\r\n";
    return resp.str() + body;
}

// =============================================================================
// Animation keyframe handlers
// =============================================================================

// POST /api/animation/save
// Body JSON: {"upload_id":"...", "fps":24, "total_frames":250, "keyframes":[...]}
// Returns: {"id":"..."}
static std::string handle_animation_save(const HttpRequest& req) {
    std::string upload_id   = json_string_value(req.body, "upload_id");
    int         fps         = json_int_value(req.body, "fps", 24);
    int         total_frames = json_int_value(req.body, "total_frames", 250);

    // Extract keyframes JSON array from body
    std::string kf_json;
    {
        std::string search = "\"keyframes\"";
        size_t pos = req.body.find(search);
        if (pos != std::string::npos) {
            pos += search.size();
            while (pos < req.body.size() && (req.body[pos] == ' ' || req.body[pos] == ':' || req.body[pos] == '\t'))
                ++pos;
            if (pos < req.body.size() && req.body[pos] == '[') {
                int depth = 0;
                size_t start = pos;
                while (pos < req.body.size()) {
                    if (req.body[pos] == '[') ++depth;
                    if (req.body[pos] == ']') { --depth; if (depth == 0) { kf_json = req.body.substr(start, pos - start + 1); break; } }
                    ++pos;
                }
            }
        }
    }
    if (kf_json.empty()) kf_json = "[]";

    if (!SupabaseClient::isConfigured()) {
        // No Supabase — generate a local id and return it
        std::string local_id = gen_id();
        return ok_json("{\"id\":\"" + local_id + "\",\"warning\":\"Supabase not configured — keyframes not persisted\"}");
    }

    std::string new_id = gen_id();
    std::ostringstream row;
    row << "{\"id\":\""           << new_id
        << "\",\"upload_id\":\""  << json_escape(upload_id)
        << "\",\"frame_data\":\""  << json_escape(kf_json)
        << "\",\"fps\":"           << fps
        << ",\"total_frames\":"    << total_frames << "}";

    std::string inserted_id = SupabaseClient::insertRow("animation_keyframes", row.str());
    if (inserted_id.empty()) inserted_id = new_id;

    return ok_json("{\"id\":\"" + json_escape(inserted_id) + "\"}");
}

// GET /api/animation/:upload_id
// Returns the keyframe row(s) for the given upload_id.
static std::string handle_animation_get(const std::string& upload_id) {
    if (!SupabaseClient::isConfigured())
        return error_json(503, "Supabase not configured");

    std::string filter = "upload_id=eq." + upload_id;
    std::string rows   = SupabaseClient::queryRows("animation_keyframes", filter);
    if (rows.empty())
        return error_json(404, "No animation data for upload_id: " + upload_id);

    return http_response(200, "OK", "application/json", rows);
}

// =============================================================================
// Shader library handlers
// =============================================================================

// GET /api/shaders  — list all shaders as JSON array
static std::string handle_shaders_list() {
    std::string body = shaders::list_json();
    return http_response(200, "OK", "application/json", body);
}

// GET /api/shaders/:name  — single shader definition (vert + frag source)
static std::string handle_shader(const std::string& name) {
    std::string body = shaders::shader_json(name);
    if (body.empty())
        return error_json(404, "Shader not found: " + name);
    return http_response(200, "OK", "application/json", body);
}

// =============================================================================
// POST /api/sim/run  — LBM flow simulation
// Body: {nx,ny,reynolds,inlet_vel,steps,obstacles:[{cx,cy,r},...]}
// =============================================================================

// Parse a JSON array of obstacle objects from the body string
static std::vector<SimObstacle> parse_obstacles(const std::string& body) {
    std::vector<SimObstacle> out;
    size_t arr = body.find("\"obstacles\"");
    if (arr == std::string::npos) return out;
    size_t ob = body.find('[', arr);
    size_t cb = body.find(']', ob);
    if (ob == std::string::npos || cb == std::string::npos) return out;
    std::string seg = body.substr(ob, cb - ob + 1);
    // Scan for {...} objects
    size_t p = 0;
    while (true) {
        size_t lo = seg.find('{', p);
        if (lo == std::string::npos) break;
        size_t hi = seg.find('}', lo);
        if (hi == std::string::npos) break;
        std::string obj = seg.substr(lo, hi - lo + 1);
        SimObstacle o;
        o.cx = (float)json_double_value(obj, "cx", 0.4);
        o.cy = (float)json_double_value(obj, "cy", 0.5);
        o.r  = (float)json_double_value(obj, "r",  0.06);
        out.push_back(o);
        p = hi + 1;
    }
    return out;
}

static std::string handle_sim_run(const HttpRequest& req) {
    SimParams sp;
    sp.nx        = json_int_value   (req.body, "nx",        200);
    sp.ny        = json_int_value   (req.body, "ny",         80);
    sp.reynolds  = (float)json_double_value(req.body, "reynolds",  150.0);
    sp.inlet_vel = (float)json_double_value(req.body, "inlet_vel",  0.1);
    sp.steps     = json_int_value   (req.body, "steps",    4000);
    sp.obstacles = parse_obstacles(req.body);

    // Clamp grid to reasonable limits
    sp.nx = std::max(40, std::min(sp.nx, 400));
    sp.ny = std::max(20, std::min(sp.ny, 200));
    sp.steps = std::max(500, std::min(sp.steps, 20000));

    SimResult res = run_lbm(sp);

    // Encode result as JSON
    // Arrays encoded compactly with 4-decimal precision
    auto arr_f = [](const std::vector<float>& v, int prec) {
        std::ostringstream os;
        os << '[';
        for (size_t i = 0; i < v.size(); ++i) {
            if (i) os << ',';
            // Fixed precision for compactness
            char buf[32];
            snprintf(buf, sizeof(buf), prec==4 ? "%.4g" : "%.3g", v[i]);
            os << buf;
        }
        os << ']';
        return os.str();
    };
    auto arr_u8 = [](const std::vector<uint8_t>& v) {
        std::ostringstream os;
        os << '[';
        for (size_t i = 0; i < v.size(); ++i) { if(i) os<<','; os<<(int)v[i]; }
        os << ']';
        return os.str();
    };

    std::ostringstream js;
    js << "{\"nx\":" << res.nx
       << ",\"ny\":" << res.ny
       << ",\"steps_done\":" << res.steps_done
       << ",\"vx\":"    << arr_f(res.vx, 4)
       << ",\"vy\":"    << arr_f(res.vy, 4)
       << ",\"speed\":" << arr_f(res.speed, 4)
       << ",\"solid\":" << arr_u8(res.solid)
       << "}";
    return ok_json(js.str());
}

// =============================================================================
// POST /api/sim/streamlines  — 3D potential-flow streamline computation
// Body: {objects:[{cx,cy,cz,sx,sy,sz},...], flow_dir:[1,0,0],
//        inlet_speed, seed_count, max_steps}
// =============================================================================

static std::vector<FlowBox> parse_flow_boxes(const std::string& body) {
    std::vector<FlowBox> out;
    size_t arr = body.find("\"objects\"");
    if (arr == std::string::npos) return out;
    size_t ob = body.find('[', arr);
    size_t cb = body.find(']', ob);
    if (ob == std::string::npos || cb == std::string::npos) return out;
    std::string seg = body.substr(ob, cb - ob + 1);
    size_t p = 0;
    while (true) {
        size_t lo = seg.find('{', p);
        if (lo == std::string::npos) break;
        size_t hi = seg.find('}', lo);
        if (hi == std::string::npos) break;
        std::string obj = seg.substr(lo, hi - lo + 1);
        FlowBox b;
        b.cx = (float)json_double_value(obj,"cx",0);
        b.cy = (float)json_double_value(obj,"cy",0);
        b.cz = (float)json_double_value(obj,"cz",0);
        b.sx = (float)json_double_value(obj,"sx",1);
        b.sy = (float)json_double_value(obj,"sy",1);
        b.sz = (float)json_double_value(obj,"sz",1);
        out.push_back(b);
        p = hi + 1;
    }
    return out;
}

static std::string handle_sim_streamlines(const HttpRequest& req) {
    FlowParams fp;
    fp.objects     = parse_flow_boxes(req.body);
    fp.inlet_speed = (float)json_double_value(req.body, "inlet_speed", 1.0);
    fp.seed_count  = std::max(10, std::min(json_int_value(req.body,"seed_count",200), 600));
    fp.max_steps   = std::max(50, std::min(json_int_value(req.body,"max_steps",400), 1000));

    // Parse flow_dir array
    {
        size_t pos = req.body.find("\"flow_dir\"");
        if (pos != std::string::npos) {
            size_t lb = req.body.find('[', pos);
            size_t rb = req.body.find(']', lb);
            if (lb != std::string::npos && rb != std::string::npos) {
                std::string seg = req.body.substr(lb+1, rb-lb-1);
                // parse 3 numbers separated by commas
                float vals[3] = {1,0,0};
                int vi = 0;
                size_t p2 = 0;
                while (vi < 3) {
                    while (p2 < seg.size() && (seg[p2]==' '||seg[p2]==',')) p2++;
                    if (p2 >= seg.size()) break;
                    vals[vi++] = std::stof(seg.substr(p2));
                    while (p2 < seg.size() && seg[p2] != ',') p2++;
                }
                memcpy(fp.flow_dir, vals, 12);
            }
        }
    }

    FlowResult res = compute_flow(fp);

    // Encode result as compact JSON
    std::ostringstream js;
    js << "{\"vmax\":" << res.vmax
       << ",\"domain\":[" << res.domain[0] << ',' << res.domain[1] << ',' << res.domain[2]
                  << ',' << res.domain[3] << ',' << res.domain[4] << ',' << res.domain[5] << ']'
       << ",\"lines\":[";

    char buf[32];
    bool first_line = true;
    for (const auto& line : res.lines) {
        if (!first_line) js << ',';
        first_line = false;
        js << "{\"p\":[";
        for (size_t i = 0; i < line.pts.size(); ++i) {
            if (i) js << ',';
            snprintf(buf, sizeof(buf), "%.4g", line.pts[i]);
            js << buf;
        }
        js << "],\"s\":[";
        for (size_t i = 0; i < line.speed.size(); ++i) {
            if (i) js << ',';
            snprintf(buf, sizeof(buf), "%.4g", line.speed[i]);
            js << buf;
        }
        js << "]}";
    }
    js << "]}";
    return ok_json(js.str());
}

// =============================================================================
// POST /api/ai/analyze  — CFD mesh AI advisor
// Body: { goal, physics, flow_dir:[x,y,z], objects:[{bbox_cx,bbox_cy,bbox_cz,
//         bbox_sx,bbox_sy,bbox_sz,tri_count,filename},...] }
// =============================================================================

struct AiObject {
    double cx, cy, cz;   // bounding box center
    double sx, sy, sz;   // bounding box size
    int    tri_count;
    std::string filename;
};

// Parse objects array from AI request body
static std::vector<AiObject> ai_parse_objects(const std::string& body) {
    std::vector<AiObject> out;
    size_t arr = body.find("\"objects\"");
    if (arr == std::string::npos) return out;
    size_t ob = body.find('[', arr);
    if (ob == std::string::npos) return out;
    size_t depth = 0, p = ob;
    size_t arr_end = ob;
    for (; p < body.size(); ++p) {
        if (body[p] == '[') ++depth;
        else if (body[p] == ']') { if (--depth == 0) { arr_end = p; break; } }
    }
    // Scan for { } objects
    p = ob + 1;
    while (p < arr_end) {
        size_t lo = body.find('{', p);
        if (lo == std::string::npos || lo >= arr_end) break;
        size_t hi = body.find('}', lo);
        if (hi == std::string::npos || hi > arr_end) break;
        std::string obj = body.substr(lo, hi - lo + 1);
        AiObject o;
        o.cx        = json_double_value(obj, "bbox_cx", 0);
        o.cy        = json_double_value(obj, "bbox_cy", 0);
        o.cz        = json_double_value(obj, "bbox_cz", 0);
        o.sx        = json_double_value(obj, "bbox_sx", 1);
        o.sy        = json_double_value(obj, "bbox_sy", 1);
        o.sz        = json_double_value(obj, "bbox_sz", 1);
        o.tri_count = json_int_value   (obj, "tri_count", 0);
        o.filename  = json_string_value(obj, "filename");
        out.push_back(o);
        p = hi + 1;
    }
    return out;
}

static std::string handle_ai_analyze(const HttpRequest& req) {
    std::string goal    = json_string_value(req.body, "goal");
    std::string physics = json_string_value(req.body, "physics");
    if (goal.empty())    goal    = "balanced";
    if (physics.empty()) physics = "aerodynamics";

    std::vector<AiObject> objs = ai_parse_objects(req.body);

    // ── Scene geometry analysis ──────────────────────────────────────────────
    double scene_min_x=1e18,scene_min_y=1e18,scene_min_z=1e18;
    double scene_max_x=-1e18,scene_max_y=-1e18,scene_max_z=-1e18;
    int total_tris = 0;

    for (auto& o : objs) {
        scene_min_x = std::min(scene_min_x, o.cx - o.sx/2);
        scene_min_y = std::min(scene_min_y, o.cy - o.sy/2);
        scene_min_z = std::min(scene_min_z, o.cz - o.sz/2);
        scene_max_x = std::max(scene_max_x, o.cx + o.sx/2);
        scene_max_y = std::max(scene_max_y, o.cy + o.sy/2);
        scene_max_z = std::max(scene_max_z, o.cz + o.sz/2);
        total_tris += o.tri_count;
    }
    if (objs.empty()) {
        scene_min_x=scene_min_y=scene_min_z=0;
        scene_max_x=scene_max_y=scene_max_z=1;
    }

    double sw = scene_max_x - scene_min_x;
    double sh = scene_max_y - scene_min_y;
    double sd = scene_max_z - scene_min_z;
    double char_len = std::max({sw, sh, sd, 1e-9});
    double min_dim  = std::max(std::min({sw, sh, sd}), 1e-9);
    double aspect   = char_len / min_dim;

    // ── Feature detection ────────────────────────────────────────────────────
    // Sharp edges: high triangle density relative to surface area heuristic
    double surf_area = 0;
    for (auto& o : objs)
        surf_area += 2.0*(o.sx*o.sy + o.sy*o.sz + o.sx*o.sz);
    double tri_density = (surf_area > 1e-9) ? (total_tris / surf_area) : 0;
    bool has_sharp   = (tri_density > 30.0) || (total_tris > 20000);

    // Symmetry: check if scene height roughly equals scene width
    double sym_ratio = sw / (sh + 1e-9);
    bool has_symm    = (sym_ratio > 0.4 && sym_ratio < 2.5) && (objs.size() == 1);

    // Wake: elongated in flow direction (aspect ratio > 1.8)
    bool has_wake    = (aspect > 1.8);

    // Inlets/outlets: pipe-like geometry (elongated + near-circular cross section)
    double cross_ratio = sw / (sh + 1e-9);
    bool has_io      = (physics == "internal_flow") ||
                       (aspect > 2.5 && cross_ratio > 0.6 && cross_ratio < 1.7);

    bool has_walls   = !objs.empty();

    // ── Compute mesh recommendations ─────────────────────────────────────────
    // Base element size as fraction of characteristic length
    double base_frac;
    if      (goal == "speed")    base_frac = 0.060;
    else if (goal == "accuracy") base_frac = 0.018;
    else                          base_frac = 0.035;

    double base_size = char_len * base_frac;

    int bl_layers;
    double bl_growth, wake_factor, surface_factor;
    if (goal == "speed") {
        bl_layers = 5;  bl_growth = 1.3;  wake_factor = 0.25; surface_factor = 0.4;
    } else if (goal == "accuracy") {
        bl_layers = 15; bl_growth = 1.15; wake_factor = 0.08; surface_factor = 0.12;
    } else {
        bl_layers = 10; bl_growth = 1.2;  wake_factor = 0.15; surface_factor = 0.20;
    }

    // Physics modifiers
    std::string turb_model;
    std::string time_stepping;
    double conv_target;
    std::string re_regime;

    if (physics == "aerodynamics") {
        turb_model    = (goal == "accuracy") ? "Spalart-Allmaras" : "k-omega SST";
        time_stepping = "steady-state";
        conv_target   = (goal == "accuracy") ? 1e-6 : 1e-4;
        re_regime     = "turbulent";
        if (goal == "accuracy") { bl_layers = std::max(bl_layers, 15); bl_growth = 1.15; }
    } else if (physics == "internal_flow") {
        turb_model    = "k-epsilon RNG";
        time_stepping = "steady-state";
        conv_target   = 1e-5;
        re_regime     = "turbulent";
        bl_layers    += 3;
    } else if (physics == "heat_transfer") {
        turb_model    = "k-omega SST";
        time_stepping = "transient";
        conv_target   = 1e-6;
        re_regime     = "laminar-turbulent";
        bl_layers     = std::max(bl_layers + 5, 15);
        bl_growth     = std::min(bl_growth, 1.15);
        surface_factor *= 0.5;
    } else if (physics == "combustion") {
        turb_model    = "k-epsilon Realizable";
        time_stepping = "transient";
        conv_target   = 1e-5;
        re_regime     = "reacting";
        base_size    *= 0.7;
        bl_layers    += 2;
    }

    double wake_size    = base_size * wake_factor;
    double surface_size = base_size * surface_factor;

    // ── Build recommendations list ───────────────────────────────────────────
    std::vector<std::string> recs, warns;

    if (physics == "aerodynamics") {
        if (has_wake) {
            recs.push_back("Extend wake refinement zone 8x chord length downstream of the trailing edge");
            recs.push_back("Use progressive coarsening in the far-field (3 refinement levels)");
        }
        if (has_symm)
            recs.push_back("Apply symmetry boundary on XZ plane to halve cell count without accuracy loss");
        if (has_sharp)
            recs.push_back("Refine leading/trailing edges with surface size reduced to " +
                           std::to_string((int)(surface_size*1000)/1.0) + " mm — critical for lift/drag accuracy");
        recs.push_back("Target y+ < 1 at all wall surfaces for wall-resolved turbulence treatment");
        recs.push_back("Domain should extend 20x chord length in all directions from the object");
    } else if (physics == "internal_flow") {
        if (has_io)
            recs.push_back("Set inlet as velocity inlet, outlet as pressure outlet with zero gauge pressure");
        recs.push_back("Use structured hex mesh in straight duct sections — reduces numerical diffusion by 40%");
        recs.push_back("Refine mesh in bends and junctions with a 3x local size reduction");
        if (has_sharp)
            recs.push_back("Apply fillets to sharp internal corners; if not possible, add local refinement");
    } else if (physics == "heat_transfer") {
        recs.push_back("Resolve thermal boundary layer: first cell height = " +
                       std::to_string((int)(surface_size*10000)/10.0) + " mm");
        recs.push_back("Enable conjugate heat transfer if solid conduction is present");
        recs.push_back("Use second-order upwind scheme for energy equation");
        if (bl_layers >= 15)
            recs.push_back("Inflation layers critical — " + std::to_string(bl_layers) +
                           " layers at growth rate " + std::to_string((int)(bl_growth*100)/100.0));
    } else if (physics == "combustion") {
        recs.push_back("Apply very fine mesh in flame zone (local size 0.5x base)");
        recs.push_back("Use species transport with finite-rate/eddy-dissipation model");
        recs.push_back("Ensure mesh resolves stoichiometric mixture fraction gradients");
        warns.push_back("Combustion simulation requires transient solver — steady-state may not converge");
    }

    // General recommendations
    if (goal == "accuracy")
        recs.push_back("Run mesh independence study: coarsen by 1.5x and verify results change < 2%");
    if (total_tris > 50000)
        warns.push_back("High source triangle count (" + std::to_string(total_tris) +
                        ") — consider decimating input STL to improve mesh quality");
    if (aspect > 5.0)
        warns.push_back("Extreme aspect ratio (" + std::to_string((int)aspect) +
                        "x) — ensure domain is long enough to capture full wake");

    // ── Build JSON response ──────────────────────────────────────────────────
    std::ostringstream o;
    o << std::fixed;

    o << "{\"features\":{"
      << "\"sharp_edges\":"    << (has_sharp ? "true" : "false") << ","
      << "\"inlets_outlets\":" << (has_io    ? "true" : "false") << ","
      << "\"walls\":"          << (has_walls ? "true" : "false") << ","
      << "\"symmetry_planes\":" << (has_symm ? "true" : "false") << ","
      << "\"wake_regions\":"   << (has_wake  ? "true" : "false")
      << "},";

    o << "\"mesh\":{"
      << "\"base_size\":"          << base_size         << ","
      << "\"surface_size\":"       << surface_size      << ","
      << "\"wake_size\":"          << wake_size         << ","
      << "\"boundary_layers\":"    << bl_layers         << ","
      << "\"bl_growth_rate\":"     << bl_growth         << ","
      << "\"char_length\":"        << char_len          << ","
      << "\"total_cells_est\":"    << (int)(std::pow(char_len/base_size,3)*0.15)
      << "},";

    o << "\"solver\":{"
      << "\"turbulence_model\":\""  << json_escape(turb_model)    << "\","
      << "\"time_stepping\":\""     << json_escape(time_stepping) << "\","
      << "\"convergence_target\":"  << conv_target                << ","
      << "\"reynolds_regime\":\""   << json_escape(re_regime)     << "\""
      << "},";

    // Recommendations array
    o << "\"recommendations\":[";
    for (size_t i = 0; i < recs.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(recs[i]) << "\"";
    }
    o << "],";

    // Warnings array
    o << "\"warnings\":[";
    for (size_t i = 0; i < warns.size(); ++i) {
        if (i) o << ",";
        o << "\"" << json_escape(warns[i]) << "\"";
    }
    o << "],";

    // Summary string
    std::string summary = "Analyzed " + std::to_string(objs.size()) + " object(s), " +
        std::to_string(total_tris) + " triangles. " +
        "Goal: " + goal + " · Physics: " + physics + ". " +
        "Recommended base size: " + std::to_string(base_size).substr(0,6) + " · " +
        std::to_string(bl_layers) + " boundary layers · " + turb_model;
    o << "\"summary\":\"" << json_escape(summary) << "\""
      << "}";

    return ok_json(o.str());
}

// =============================================================================
// POST /api/ai/chat — CFD knowledge-base chat engine
// Body: { message, goal, physics }
// =============================================================================
static std::string handle_ai_chat(const HttpRequest& req) {
    std::string msg     = json_string_value(req.body, "message");
    std::string goal    = json_string_value(req.body, "goal");
    std::string physics = json_string_value(req.body, "physics");
    if (msg.empty()) return error_json(400, "message required");

    // Lowercase for matching
    std::string ml = msg;
    std::transform(ml.begin(), ml.end(), ml.begin(), ::tolower);

    std::string reply;
    std::vector<std::string> sugg;

    auto contains = [&](const std::string& s){ return ml.find(s) != std::string::npos; };

    if (contains("y+") || contains("y plus") || contains("yplus") || contains("wall distance")) {
        reply = "y+ is the dimensionless wall distance critical for turbulence accuracy.\n\n"
                "• Wall-resolved (k-ω SST): y+ < 1 — resolves the viscous sublayer\n"
                "• Wall functions: 30 < y+ < 300 — uses log-law approximation\n\n"
                "First cell height h ≈ (y+ × ν) / u_τ\n"
                "where u_τ ≈ 0.05 × U_∞ for external flow.\n\n"
                "For " + (physics.empty() ? "aerodynamics" : physics) + " with " +
                (goal == "accuracy" ? "accuracy goal, target y+ < 1." :
                 goal == "speed"    ? "speed goal, y+ 30–100 with wall functions is acceptable." :
                                     "balanced goal, aim for y+ < 5 with enhanced wall treatment.");
        sugg = {"How to set first cell height in mesh?","Wall functions vs wall-resolved — which to choose?","k-omega SST y+ requirements"};

    } else if (contains("turbulence") || contains("k-omega") || contains("k-epsilon") || contains("spalart") || contains("rans")) {
        if      (physics == "aerodynamics")   reply = "For aerodynamics: k-ω SST is the gold standard. It switches between k-ε (freestream) and k-ω (near wall), handling adverse pressure gradients and trailing-edge separation well.\n\nSpalart-Allmaras is faster and works well for attached flows (wings at low AoA).";
        else if (physics == "internal_flow")  reply = "For internal flow: k-ε Realizable handles recirculation, jets, and separating flow better than standard k-ε. k-ε RNG is good when swirl is present (cyclones, curved ducts).";
        else if (physics == "heat_transfer")  reply = "For conjugate heat transfer: k-ω SST with Enhanced Wall Treatment. It accurately resolves the thermal boundary layer when y+ < 1. The Prandtl turbulent number Pr_t ≈ 0.85 is used for the energy equation.";
        else if (physics == "combustion")     reply = "For combustion: k-ε Realizable with Eddy Dissipation Concept (EDC) for detailed chemistry, or simple Eddy Dissipation for fast reactions. Consider enabling Radiation (P1 or Discrete Ordinates) for high-temperature flames.";
        else                                  reply = "Turbulence model guide:\n• k-ω SST — best all-round for external aero\n• k-ε Realizable — internal flows, jets\n• k-ε RNG — swirling/rotating flows\n• Spalart-Allmaras — simple aero, fast\n• RSM — highly anisotropic flows (expensive)";
        sugg = {"What y+ for this model?","Inlet turbulence intensity values","When does turbulence model choice matter most?"};

    } else if (contains("mesh") && (contains("size") || contains("refine") || contains("coarse") || contains("fine"))) {
        reply = "Mesh sizing rules of thumb:\n\n"
                "Base size (% of char. length):\n"
                "• Speed goal: 5–7%\n"
                "• Balanced: 3–4%\n"
                "• Accuracy: 1–2%\n\n"
                "Refinement zones:\n"
                "• Surface: 20–40% of base size\n"
                "• Wake: 8–15% of base, extend 5–8× downstream\n"
                "• Boundary layers: 10–15 layers, growth rate 1.15–1.2\n\n"
                "Click ⚡ Analyze for scene-specific values.";
        sugg = {"How many boundary layers?","Wake zone dimensions","Mesh independence study"};

    } else if (contains("boundary condition") || contains("inlet") || contains("outlet") || contains("bc ")) {
        reply = "Standard boundary conditions:\n\n"
                "• Velocity inlet — specify U, k, ε (or ω), temperature\n"
                "• Pressure outlet — gauge pressure = 0 Pa (recommended)\n"
                "• No-slip wall — default; add heat flux for thermal problems\n"
                "• Symmetry — zero normal flux; valid only when flow is symmetric\n"
                "• Periodic — for repeating geometries (blade passages, fins)\n"
                "• Farfield — for external compressible aerodynamics\n\n"
                "Domain size: place inlet/outlet at least 10× L from the object.";
        sugg = {"Turbulence intensity at inlet","Domain sizing for external flow","When to use symmetry?"};

    } else if (contains("converg") || contains("residual") || contains("diverge") || contains("not converging")) {
        reply = "Convergence troubleshooting:\n\n"
                "1. Check mesh quality — skewness < 0.85, aspect ratio < 100\n"
                "2. Reduce under-relaxation: pressure 0.3, momentum 0.5\n"
                "3. Start with 1st-order schemes, switch to 2nd after ~200 iter.\n"
                "4. Initialize from freestream, not zero\n"
                "5. Monitor Cl/Cd alongside residuals — they must stabilize too\n"
                "6. If oscillating: flow may be unsteady — switch to transient\n\n"
                "Target: continuity < 1e-4 for engineering, < 1e-6 for research.";
        sugg = {"How to check mesh quality?","Switching to transient solver","Under-relaxation factor guide"};

    } else if (contains("lift") || contains("drag") || contains(" cd ") || contains(" cl ") || contains("coefficient")) {
        reply = "Aerodynamic coefficients:\n\n"
                "Cd = F_drag / (½ ρ U² A_ref)\n"
                "Cl = F_lift / (½ ρ U² A_ref)\n\n"
                "Best practices:\n"
                "• Use 2nd-order upwind for momentum — critical for force accuracy\n"
                "• Pressure + viscous contributions reported separately\n"
                "• A_ref: frontal area (bluff bodies) or planform area (wings)\n"
                "• Run 500–2000 iterations and average last 200 if oscillating\n"
                "• Mesh sensitivity: Cd should change < 1% on refinement";
        sugg = {"Reference area setup","Pressure vs viscous drag breakdown","How to reduce drag numerically?"};

    } else if (contains("heat") || contains("thermal") || contains("temperature") || contains("nusselt")) {
        reply = "Heat transfer setup:\n\n"
                "• Enable energy equation in solver\n"
                "• Wall BC: constant T (isothermal) or constant q'' (heat flux)\n"
                "• Nu = h·L/k — post-process from wall heat flux\n"
                "• Resolve thermal BL: y+ < 1 mandatory for accuracy\n"
                "• Pr (water) ≈ 7, Pr (air) ≈ 0.71, Pr (oil) >> 1\n"
                "• Thinner thermal BL for high Pr fluids — needs finer wall mesh";
        sugg = {"Conjugate heat transfer setup","Convective heat transfer coefficient","Thermal boundary layer thickness"};

    } else if (contains("domain") || contains("far field") || contains("farfield") || contains("domain size")) {
        std::string sz = (physics == "aerodynamics") ? "20–30× chord" :
                         (physics == "internal_flow") ? "10–20× diameter for inlet/outlet extensions" :
                                                       "15× characteristic length";
        reply = "Domain sizing for " + (physics.empty() ? "external flow" : physics) + ":\n\n"
                "• Recommended extent: " + sz + "\n"
                "• Upstream: 5–8× L (less critical)\n"
                "• Downstream: 15–20× L (captures wake)\n"
                "• Lateral: 10× L (minimizes blockage effects)\n\n"
                "Blockage ratio = A_object / A_domain cross-section < 5% for valid results.\n"
                "Higher blockage artificially increases drag.";
        sugg = {"What is blockage ratio?","Inlet placement","Outlet backflow issues"};

    } else if (contains("what") && (contains("do") || contains("how") || contains("should") || contains("recommend"))) {
        reply = "Here's how to get started with your CFD analysis:\n\n"
                "1. Select your Goal (Speed / Balanced / Accuracy) above\n"
                "2. Select the Physics Type matching your application\n"
                "3. Click ⚡ Detect to auto-identify geometric features\n"
                "4. Click Analyze — the C++ engine will compute mesh settings,\n"
                "   boundary layer parameters, and solver recommendations\n\n"
                "Then ask me any specific question about your setup!";
        sugg = {"What turbulence model for aerodynamics?","How many boundary layers?","Explain y+ to me"};

    } else {
        reply = "I'm the Discreetize CFD advisor powered by the C++ analysis engine. I can answer questions about:\n\n"
                "• Turbulence models (k-ω SST, k-ε, SA, RSM)\n"
                "• Mesh settings (sizing, boundary layers, wake zones)\n"
                "• Boundary conditions (inlet, outlet, symmetry, walls)\n"
                "• Convergence & solver settings\n"
                "• Physics (lift/drag, heat transfer, combustion)\n"
                "• y+ and wall treatment\n\n"
                "Try the Analyze button for scene-specific recommendations!";
        sugg = {"Best turbulence model for my physics?","How to set up boundary layers?","Why isn't my simulation converging?"};
    }

    std::ostringstream o;
    o << "{\"reply\":\"" << json_escape(reply) << "\",\"suggestions\":[";
    for (size_t i = 0; i < sugg.size(); ++i) { if (i) o << ","; o << "\"" << json_escape(sugg[i]) << "\""; }
    o << "]}";
    return ok_json(o.str());
}

// =============================================================================
// Request dispatcher
// =============================================================================

static std::string dispatch(const HttpRequest& req) {
    const std::string& method = req.method;
    std::string        route  = strip_query(req.path);

    // Health check
    if (method == "GET" && (route == "/health" || route == "/api/health"))
        return handle_health();

    // Upload
    if (method == "POST" && route == "/api/upload")
        return handle_upload(req);

    // Generate mesh
    if (method == "POST" && route == "/api/generate-mesh")
        return handle_generate_mesh(req);

    // List all uploads: GET /api/uploads
    if (method == "GET" && route == "/api/uploads")
        return handle_list_uploads();

    // Upload sub-routes: GET /api/upload/<id>/faces|stl
    //                    DELETE /api/upload/<id>
    if (starts_with(route, "/api/upload/")) {
        std::string rest = path_after(route, "/api/upload/");
        if (method == "DELETE") {
            // rest is just the upload ID (no trailing slash expected)
            std::string upload_id = rest;
            if (upload_id.empty())
                return error_json(400, "Missing upload_id");
            return handle_delete_upload(upload_id);
        }
        if (method == "GET") {
            size_t slash = rest.find('/');
            if (slash != std::string::npos) {
                std::string upload_id = rest.substr(0, slash);
                std::string endpoint  = rest.substr(slash + 1);
                if (endpoint == "faces") return handle_upload_faces(upload_id);
                if (endpoint == "stl")   return handle_upload_stl(upload_id);
            }
            return error_json(404, "Unknown upload endpoint");
        }
    }

    // Mesh sub-routes: GET /api/mesh/<id>/{status,data,stl}
    if (method == "GET" && starts_with(route, "/api/mesh/")) {
        std::string rest = path_after(route, "/api/mesh/");
        size_t slash = rest.find('/');
        if (slash != std::string::npos) {
            std::string mesh_id  = rest.substr(0, slash);
            std::string endpoint = rest.substr(slash + 1);
            if (endpoint == "status") return handle_mesh_status(mesh_id);
            if (endpoint == "data")   return handle_mesh_data(mesh_id);
            if (endpoint == "stl")    return handle_mesh_stl(mesh_id);
        }
        return error_json(404, "Unknown mesh endpoint");
    }

    // Primitives: POST /api/primitives/box|sphere
    if (method == "POST" && starts_with(route, "/api/primitives/")) {
        std::string type = path_after(route, "/api/primitives/");
        if (!type.empty()) return handle_primitive(type, req);
        return error_json(400, "Missing primitive type");
    }

    // Boolean subtract: POST /api/boolean/subtract (box-ellipsoid, legacy)
    if (method == "POST" && route == "/api/boolean/subtract")
        return handle_boolean_subtract(req);

    // ── v1 API ───────────────────────────────────────────────────────────────
    // Mesh-mesh boolean subtraction (async job)
    if (method == "POST" && route == "/api/v1/boolean/subtract-meshes")
        return handle_v1_boolean_subtract(req);

    if (starts_with(route, "/api/v1/boolean/jobs/")) {
        std::string rest = path_after(route, "/api/v1/boolean/jobs/");
        size_t slash = rest.find('/');
        std::string jid  = (slash==std::string::npos) ? rest : rest.substr(0,slash);
        std::string ep   = (slash==std::string::npos) ? "" : rest.substr(slash+1);
        if(jid.empty()) return error_json(400,"missing job id");
        if(method=="GET" && ep.empty())          return handle_v1_boolean_status(jid);
        if(method=="GET" && ep=="result.stl")    return handle_v1_boolean_stl(jid);
        if(method=="GET" && ep=="result.vtu")    return handle_v1_boolean_vtu(jid);
        return error_json(404,"unknown boolean job endpoint");
    }

    // Bifilar solid (async job)
    if (method == "POST" && route == "/api/v1/bifilar/generate")
        return handle_v1_bifilar_generate(req);

    if (starts_with(route, "/api/v1/bifilar/jobs/")) {
        std::string rest = path_after(route, "/api/v1/bifilar/jobs/");
        size_t slash = rest.find('/');
        std::string jid  = (slash==std::string::npos) ? rest : rest.substr(0,slash);
        std::string ep   = (slash==std::string::npos) ? "" : rest.substr(slash+1);
        if(jid.empty()) return error_json(400,"missing job id");
        if(method=="GET" && ep.empty())          return handle_v1_bifilar_status(jid);
        if(method=="GET" && ep=="result.stl")    return handle_v1_bifilar_stl(jid);
        if(method=="GET" && ep=="result.vtu")    return handle_v1_bifilar_vtu(jid);
        return error_json(404,"unknown bifilar job endpoint");
    }

    if (method == "POST" && route == "/api/v1/bifilar/generate-visual")
        return handle_v1_bifilar_visual(req);

    if (starts_with(route, "/api/v1/bifilar/visual-jobs/")) {
        std::string rest = path_after(route, "/api/v1/bifilar/visual-jobs/");
        size_t slash = rest.find('/');
        std::string jid = (slash==std::string::npos) ? rest : rest.substr(0,slash);
        std::string ep  = (slash==std::string::npos) ? "" : rest.substr(slash+1);
        if(jid.empty()) return error_json(400,"missing job id");
        if(method=="GET" && ep.empty())                       return handle_v1_bifilar_visual_status(jid);
        if(method=="GET" && (ep=="inner1.stl"||ep=="inner2.stl")) return handle_v1_bifilar_visual_stl(jid,ep);
        return error_json(404,"unknown visual-job endpoint");
    }

    // Task pipeline triggers
    if (method == "GET"  && route == "/api/v1/task1/models") return handle_v1_task1_models();
    if (method == "POST" && route == "/api/v1/task1/run")    return handle_v1_task1_run(req);
    if (method == "POST" && route == "/api/v1/task2/run")    return handle_v1_task2_run(req);

    // Heal: POST /api/heal
    if (method == "POST" && route == "/api/heal")
        return handle_heal(req);

    // Slice: POST /api/slice
    if (method == "POST" && route == "/api/slice")
        return handle_slice(req);

    // Curvature: GET /api/curvature/<id>
    if (method == "GET" && starts_with(route, "/api/curvature/")) {
        std::string upload_id = path_after(route, "/api/curvature/");
        if (!upload_id.empty()) return handle_curvature(upload_id);
        return error_json(404, "Missing upload_id");
    }

    // Export: GET /api/export/<id>/<format>
    if (method == "GET" && starts_with(route, "/api/export/")) {
        std::string rest  = path_after(route, "/api/export/");
        size_t slash = rest.find('/');
        if (slash != std::string::npos) {
            std::string uid = rest.substr(0, slash);
            std::string fmt = rest.substr(slash + 1);
            if (!uid.empty() && !fmt.empty()) return handle_export(uid, fmt);
        }
        return error_json(400, "Usage: GET /api/export/<id>/<format>");
    }

    // Render
    if (method == "GET"  && route == "/api/render/defaults")
        return handle_render_defaults();
    if (method == "POST" && route == "/api/render")
        return handle_render(req);
    if (method == "POST" && route == "/api/render/frame")
        return handle_render_frame(req);

    // Shader library
    if (method == "GET" && route == "/api/shaders")
        return handle_shaders_list();
    // Node graph → GLSL compiler (must come before the generic /api/shaders/:name handler)
    if (method == "POST" && route == "/api/shaders/from_nodes")
        return handle_compile_nodes(req);
    if (method == "GET" && starts_with(route, "/api/shaders/")) {
        std::string name = path_after(route, "/api/shaders/");
        return handle_shader(name);
    }

    // Animation keyframes
    if (method == "POST" && route == "/api/animation/save")
        return handle_animation_save(req);
    if (method == "GET" && starts_with(route, "/api/animation/")) {
        std::string uid = path_after(route, "/api/animation/");
        if (!uid.empty()) return handle_animation_get(uid);
        return error_json(400, "Missing upload_id");
    }

    // CFD flow simulation (2D LBM)
    if (method == "POST" && route == "/api/sim/run")
        return handle_sim_run(req);

    // 3D potential-flow streamlines
    if (method == "POST" && route == "/api/sim/streamlines")
        return handle_sim_streamlines(req);

    if (method == "POST" && route == "/api/ai/analyze")
        return handle_ai_analyze(req);

    if (method == "POST" && route == "/api/ai/chat")
        return handle_ai_chat(req);

    // CORS preflight
    if (method == "OPTIONS")
        return http_response(204, "No Content", "text/plain", "");

    // Public config (Supabase URL + anon key for the frontend)
    if (method == "GET" && route == "/api/config") {
        const char* su = std::getenv("SUPABASE_URL");
        const char* ak = std::getenv("SUPABASE_ANON_KEY");
        std::string body = "{\"supabaseUrl\":\"";
        body += su ? json_escape(su) : "";
        body += "\",\"supabaseAnonKey\":\"";
        body += ak ? json_escape(ak) : "";
        body += "\"}";
        return ok_json(body);
    }

    // Static files / SPA fallback
    if (method == "GET")
        return handle_static(route);

    return error_json(405, "Method Not Allowed");
}

// =============================================================================
// Hot-reload file watcher
// =============================================================================

// Runs in a detached background thread.  Polls stat() on every file under
// WEBUI_ROOT every 300 ms; bumps g_reload_gen on any mtime change so that
// open SSE connections see the new generation and push a reload event.
static void hotreload_watcher() {
    std::string path = WEBUI_ROOT + "/index.html";
    struct stat st{};
    time_t last = (::stat(path.c_str(), &st) == 0) ? st.st_mtime : 0;

    while (true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        struct stat cur{};
        if (::stat(path.c_str(), &cur) == 0 && cur.st_mtime != last) {
            last = cur.st_mtime;
            ++g_reload_gen;
            std::cout << "[hotreload] index.html changed — reloading browser\n";
        }
    }
}

// =============================================================================
// Per-connection handler
// =============================================================================

static void handle_client(int client_fd) {
    // Apply a read timeout so idle connections don't block forever.
    // 300 seconds is enough for a 50 MB upload even on a slow connection
    // (e.g. 50 MB / 300 s ≈ 1.4 Mbit/s minimum throughput).
    struct timeval tv;
    tv.tv_sec  = 300;
    tv.tv_usec = 0;
    setsockopt(client_fd, SOL_SOCKET, SO_RCVTIMEO,
               reinterpret_cast<const void*>(&tv), sizeof(tv));

    std::string raw = recv_all(client_fd);
    std::string response;

    if (raw.empty()) {
        response = error_json(400, "Empty request");
    } else {
        HttpRequest req = parse_request(raw);
        if (!req.valid) {
            response = error_json(400, "Malformed HTTP request");
        } else {

            // ── Hot-reload SSE endpoint ───────────────────────────────────────
            // Keeps the socket open; bypasses the normal string-response path.
            if (req.method == "GET" && strip_query(req.path) == "/api/livereload") {
                const std::string hdr =
                    "HTTP/1.1 200 OK\r\n"
                    "Content-Type: text/event-stream\r\n"
                    "Cache-Control: no-cache\r\n"
                    "Connection: keep-alive\r\n"
                    "Access-Control-Allow-Origin: *\r\n"
                    "\r\n"
                    "data: connected\n\n";
                send(client_fd, hdr.c_str(), hdr.size(), MSG_NOSIGNAL);

                uint64_t seen = g_reload_gen.load();
                while (true) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(250));
                    uint64_t cur = g_reload_gen.load();
                    if (cur != seen) {
                        seen = cur;
                        const char* msg = "data: reload\n\n";
                        if (::send(client_fd, msg, 14, MSG_NOSIGNAL) <= 0) break;
                    }
                }
                close(client_fd);
                return;
            }
            // ─────────────────────────────────────────────────────────────────

            try {
                response = dispatch(req);
            } catch (const std::exception& ex) {
                response = error_json(500, std::string("Internal server error: ") + ex.what());
            } catch (...) {
                response = error_json(500, "Unknown internal server error");
            }
        }
    }

    // Send full response
    const char* ptr = response.data();
    size_t      rem = response.size();
    while (rem > 0) {
        ssize_t n = send(client_fd, ptr, rem, 0);
        if (n <= 0) break;
        ptr += n;
        rem -= static_cast<size_t>(n);
    }

    close(client_fd);
}

// =============================================================================
// main
// =============================================================================

int main(int argc, char** argv) {
    // Resolve webui path relative to the directory containing the executable
    {
        std::string exe = argv[0];
        auto slash = exe.rfind('/');
        std::string exe_dir = (slash == std::string::npos) ? "." : exe.substr(0, slash);
        WEBUI_ROOT = exe_dir + "/webui";
        // If that doesn't exist, fall back to ./webui (common when run from project root)
        struct stat st{};
        if (stat(WEBUI_ROOT.c_str(), &st) != 0)
            WEBUI_ROOT = "webui";
    }

    SupabaseClient::init();

    // Start hot-reload file watcher in background
    std::thread(hotreload_watcher).detach();
    std::cout << "Hot-reload watching: " << WEBUI_ROOT << "\n";

    int port = 8000;
    if (const char* ep = std::getenv("PORT")) {
        try { port = std::stoi(ep); } catch (...) {}
    }
    if (argc >= 2) {
        try { port = std::stoi(argv[1]); } catch (...) {}
    }

    // Create TCP/IPv4 socket
    int server_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (server_fd < 0) {
        std::cerr << "socket() failed: " << strerror(errno) << "\n";
        return 1;
    }

    // Reuse address/port so we can restart quickly
    int opt = 1;
    setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR,
               reinterpret_cast<const void*>(&opt), sizeof(opt));
#ifdef SO_REUSEPORT
    setsockopt(server_fd, SOL_SOCKET, SO_REUSEPORT,
               reinterpret_cast<const void*>(&opt), sizeof(opt));
#endif

    sockaddr_in addr{};
    addr.sin_family      = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port        = htons(static_cast<uint16_t>(port));

    if (bind(server_fd, reinterpret_cast<struct sockaddr*>(&addr), sizeof(addr)) < 0) {
        std::cerr << "bind() failed on port " << port << ": " << strerror(errno) << "\n";
        close(server_fd);
        return 1;
    }

    if (listen(server_fd, 128) < 0) {
        std::cerr << "listen() failed: " << strerror(errno) << "\n";
        close(server_fd);
        return 1;
    }

    std::cout << "CFD Mesh Server running at http://localhost:" << port << "\n";
    std::cout << "WebUI root: " << WEBUI_ROOT << "\n";

    while (true) {
        sockaddr_in  client_addr{};
        socklen_t    client_len = sizeof(client_addr);
        int client_fd = accept(server_fd,
                               reinterpret_cast<struct sockaddr*>(&client_addr),
                               &client_len);
        if (client_fd < 0) {
            if (errno == EINTR) continue;
            std::cerr << "accept() failed: " << strerror(errno) << "\n";
            continue;
        }

        std::thread([client_fd]() {
            handle_client(client_fd);
        }).detach();
    }

    close(server_fd);
    return 0;
}
