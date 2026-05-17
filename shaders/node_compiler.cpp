// =============================================================================
// node_compiler.cpp — Node graph → GLSL compiler
// =============================================================================
//
// Supabase SQL schema required (run in Supabase SQL editor):
//
// CREATE TABLE uploads (
//     id TEXT PRIMARY KEY,
//     filename TEXT,
//     tri_count INT,
//     created_at TIMESTAMPTZ DEFAULT NOW()
// );
//
// CREATE TABLE node_shaders (
//     id TEXT PRIMARY KEY DEFAULT gen_random_uuid(),
//     upload_id TEXT,
//     graph_json TEXT,
//     vertex_glsl TEXT,
//     fragment_glsl TEXT,
//     uniforms_json TEXT,
//     created_at TIMESTAMPTZ DEFAULT NOW()
// );
//
// CREATE TABLE animation_keyframes (
//     id TEXT PRIMARY KEY DEFAULT gen_random_uuid(),
//     upload_id TEXT,
//     frame_data TEXT,
//     fps INT DEFAULT 24,
//     total_frames INT DEFAULT 250,
//     created_at TIMESTAMPTZ DEFAULT NOW()
// );
//
// Storage buckets needed: "stl-files", "renders"
// =============================================================================

#include "node_compiler.hpp"
#include <array>
#include <cstdlib>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace node_compiler {

// =============================================================================
// Minimal JSON helpers — no external library
// =============================================================================

// Extract a string value for a key in flat JSON: "key":"value"
static std::string jstr(const std::string& json, const std::string& key,
                         const std::string& def = "") {
    std::string search = "\"" + key + "\"";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return def;
    pos += search.size();
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == ':' || json[pos] == '\t'))
        ++pos;
    if (pos >= json.size() || json[pos] != '"') return def;
    ++pos;
    std::string val;
    while (pos < json.size() && json[pos] != '"') {
        if (json[pos] == '\\' && pos + 1 < json.size()) {
            char esc = json[++pos];
            if      (esc == '"')  val += '"';
            else if (esc == '\\') val += '\\';
            else if (esc == 'n')  val += '\n';
            else if (esc == 'r')  val += '\r';
            else if (esc == 't')  val += '\t';
            else                  val += esc;
            ++pos;
        } else {
            val += json[pos++];
        }
    }
    return val;
}

// Extract a numeric value for a key: "key":number
static double jnum(const std::string& json, const std::string& key, double def = 0.0) {
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

// Extract a [r,g,b] array for a key.
static std::array<double, 3> json_vec3_val(const std::string& json, const std::string& key,
                                            std::array<double, 3> def) {
    std::string search = "\"" + key + "\"";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return def;
    pos += search.size();
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == ':' || json[pos] == '\t'))
        ++pos;
    if (pos >= json.size() || json[pos] != '[') return def;
    ++pos;
    std::array<double, 3> out{};
    for (int i = 0; i < 3; ++i) {
        while (pos < json.size() && (json[pos] == ' ' || json[pos] == ',' ||
               json[pos] == '\t' || json[pos] == '\n'))
            ++pos;
        if (pos >= json.size()) return def;
        try {
            size_t consumed = 0;
            out[i] = std::stod(json.substr(pos), &consumed);
            if (consumed == 0) return def;
            pos += consumed;
        } catch (...) { return def; }
    }
    return out;
}

// Return the raw JSON object string for a given key: "key":{...} or "key":[...]
// Returns the substring from '{' or '[' to its matching closing bracket.
static std::string jobj(const std::string& json, const std::string& key) {
    std::string search = "\"" + key + "\"";
    size_t pos = json.find(search);
    if (pos == std::string::npos) return "";
    pos += search.size();
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == ':' || json[pos] == '\t'))
        ++pos;
    if (pos >= json.size()) return "";
    char open = json[pos];
    char close = (open == '{') ? '}' : (open == '[') ? ']' : '\0';
    if (close == '\0') return "";
    int depth = 0;
    size_t start = pos;
    while (pos < json.size()) {
        if (json[pos] == open)  ++depth;
        if (json[pos] == close) { --depth; if (depth == 0) return json.substr(start, pos - start + 1); }
        ++pos;
    }
    return "";
}

// Parse array of JSON objects: finds all {...} blocks at top level within an array string.
static std::vector<std::string> parse_array_objects(const std::string& arr) {
    std::vector<std::string> result;
    size_t pos = 0;
    while (pos < arr.size()) {
        size_t start = arr.find('{', pos);
        if (start == std::string::npos) break;
        int depth = 0;
        size_t i = start;
        while (i < arr.size()) {
            if (arr[i] == '{') ++depth;
            if (arr[i] == '}') { --depth; if (depth == 0) break; }
            // Handle strings to avoid counting braces inside them
            if (arr[i] == '"') {
                ++i;
                while (i < arr.size() && arr[i] != '"') {
                    if (arr[i] == '\\') ++i;
                    ++i;
                }
            }
            ++i;
        }
        if (i < arr.size()) {
            result.push_back(arr.substr(start, i - start + 1));
        }
        pos = i + 1;
    }
    return result;
}

// =============================================================================
// GLSL source constants
// =============================================================================

static const std::string VERTEX_GLSL = R"(varying vec3 vNormal;
varying vec3 vWorldPos;
varying vec3 vViewDir;
void main(){
  vec4 wp=modelMatrix*vec4(position,1.0);
  vWorldPos=wp.xyz;
  vNormal=normalize(normalMatrix*normal);
  vViewDir=normalize(cameraPosition-wp.xyz);
  gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.0);
})";

// Full GGX PBR fragment shader with two lights, tone mapping and gamma correction.
static const std::string PBR_FRAG_TEMPLATE = R"(uniform vec3  baseColor;
uniform float roughness;
uniform float metalness;
uniform float alpha;
uniform vec3  emissionColor;
uniform float emissionStrength;

varying vec3 vNormal;
varying vec3 vWorldPos;
varying vec3 vViewDir;

const float PI = 3.14159265359;

// GGX / Trowbridge-Reitz normal distribution
float D_GGX(float NdotH, float a) {
    float a2 = a * a;
    float f  = (NdotH * a2 - NdotH) * NdotH + 1.0;
    return a2 / (PI * f * f);
}

// Smith-GGX geometry term (Schlick approximation)
float G_Smith_Schlick(float NdotV, float NdotL, float a) {
    float k   = (a + 1.0);
    k = k * k / 8.0;
    float gv  = NdotV / (NdotV * (1.0 - k) + k);
    float gl  = NdotL / (NdotL * (1.0 - k) + k);
    return gv * gl;
}

// Schlick Fresnel
vec3 F_Schlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

// ACES filmic tonemapping
vec3 aces(vec3 x) {
    const float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
    return clamp((x*(a*x+b))/(x*(c*x+d)+e), 0.0, 1.0);
}

vec3 pbr_shade(vec3 N, vec3 V, vec3 L, vec3 Lcolor,
               vec3 albedo, float rough, float metal, vec3 F0) {
    vec3  H     = normalize(V + L);
    float NdotL = max(dot(N, L), 0.0);
    float NdotV = max(dot(N, V), 1e-4);
    float NdotH = max(dot(N, H), 0.0);
    float VdotH = max(dot(V, H), 0.0);

    float a   = rough * rough;
    float D   = D_GGX(NdotH, a);
    float G   = G_Smith_Schlick(NdotV, NdotL, rough);
    vec3  F   = F_Schlick(VdotH, F0);

    vec3 kd   = (1.0 - F) * (1.0 - metal);
    vec3 spec = (D * G * F) / max(4.0 * NdotV * NdotL, 1e-4);
    vec3 diff = kd * albedo / PI;

    return (diff + spec) * Lcolor * NdotL;
}

void main() {
    vec3 N = normalize(vNormal);
    vec3 V = normalize(vViewDir);

    // Two-point light setup: key light + fill light
    vec3 L0 = normalize(vec3(2.0, 4.0, 3.0));
    vec3 L1 = normalize(vec3(-3.0, 1.0, -2.0));
    vec3 Lc0 = vec3(3.0, 2.8, 2.5);
    vec3 Lc1 = vec3(0.4, 0.45, 0.6);

    vec3 albedo = baseColor;
    float rough = clamp(roughness, 0.04, 1.0);
    float metal = clamp(metalness, 0.0,  1.0);
    vec3 F0 = mix(vec3(0.04), albedo, metal);

    vec3 Lo = pbr_shade(N, V, L0, Lc0, albedo, rough, metal, F0)
            + pbr_shade(N, V, L1, Lc1, albedo, rough, metal, F0);

    // Ambient: image-based approximation
    vec3 amb = vec3(0.03) * albedo;
    // Fresnel-weighted env reflection for metallic surfaces
    vec3 F_amb = F_Schlick(max(dot(N, V), 0.0), F0);
    vec3 env   = F_amb * (1.0 - rough) * 0.15 * albedo;

    vec3 color = Lo + amb + env + emissionColor * emissionStrength;
    color = aces(color);
    // Gamma correction
    color = pow(clamp(color, 0.0, 1.0), vec3(1.0 / 2.2));
    gl_FragColor = vec4(color, alpha);
})";

// =============================================================================
// Node data
// =============================================================================

struct NodeInfo {
    int         id   = -1;
    std::string type;
    std::string props_json;  // raw JSON of "props" object
};

struct Connection {
    int from_id   = -1;
    int from_port = 0;
    int to_id     = -1;
    int to_port   = 0;
};

struct Graph {
    std::map<int, NodeInfo> nodes;       // id → node
    std::vector<Connection> connections;
};

// =============================================================================
// Graph parsing
// =============================================================================

static Graph parse_graph(const std::string& json) {
    Graph g;

    // Parse nodes array
    std::string nodes_arr = jobj(json, "nodes");
    if (!nodes_arr.empty()) {
        auto node_objs = parse_array_objects(nodes_arr);
        for (const auto& obj : node_objs) {
            NodeInfo ni;
            ni.id   = static_cast<int>(jnum(obj, "id", -1));
            ni.type = jstr(obj, "type");
            ni.props_json = jobj(obj, "props");
            if (ni.id >= 0 && !ni.type.empty())
                g.nodes[ni.id] = ni;
        }
    }

    // Parse connections array
    std::string conns_arr = jobj(json, "connections");
    if (!conns_arr.empty()) {
        auto conn_objs = parse_array_objects(conns_arr);
        for (const auto& obj : conn_objs) {
            Connection c;
            c.from_id   = static_cast<int>(jnum(obj, "from",     -1));
            c.from_port = static_cast<int>(jnum(obj, "fromPort",  0));
            c.to_id     = static_cast<int>(jnum(obj, "to",       -1));
            c.to_port   = static_cast<int>(jnum(obj, "toPort",    0));
            if (c.from_id >= 0 && c.to_id >= 0)
                g.connections.push_back(c);
        }
    }

    return g;
}

// Find the node connected to toPort of toNode. Returns nullptr if none.
static const NodeInfo* find_input(const Graph& g, int to_id, int to_port) {
    for (const auto& c : g.connections) {
        if (c.to_id == to_id && c.to_port == to_port) {
            auto it = g.nodes.find(c.from_id);
            if (it != g.nodes.end()) return &it->second;
        }
    }
    return nullptr;
}

// =============================================================================
// Per-node shader generation
// =============================================================================

// Helper: format a double compactly (e.g. 0.5 → "0.5", not "0.500000")
static std::string fd(double v) {
    std::ostringstream ss;
    ss << v;
    std::string s = ss.str();
    // Make sure it has a decimal point so GLSL doesn't treat it as int
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos)
        s += ".0";
    return s;
}

// Helper: format a vec3 literal inline
static std::string fv3(std::array<double,3> c) {
    return "vec3(" + fd(c[0]) + "," + fd(c[1]) + "," + fd(c[2]) + ")";
}

struct ShaderResult {
    std::string frag;         // fragment GLSL body
    std::string uniforms_decl;// uniform declarations (inserted at top of frag)
    std::string uniforms_json;// JSON object of uniform values
};

// Forward declaration
static ShaderResult compile_node(const Graph& g, const NodeInfo& node,
                                  const std::string& prefix);

// ---------- Principled BSDF node -------------------------------------------

static ShaderResult compile_principled(const Graph& g, const NodeInfo& node,
                                        const std::string& pfx) {
    const std::string& props = node.props_json;

    // Defaults from props
    std::array<double,3> baseColor   = json_vec3_val(props, "baseColor",    {0.8, 0.8, 0.8});
    double roughness                 = jnum(props, "roughness",    0.5);
    double metalness                 = jnum(props, "metalness",    0.0);
    double alpha_val                 = jnum(props, "alpha",        1.0);
    std::array<double,3> emColor     = json_vec3_val(props, "emissionColor",{0.0,0.0,0.0});
    double emStrength                = jnum(props, "emissionStrength", 0.0);

    // Check if a node is overriding baseColor (port 0), roughness (port 1), metalness (port 2)
    const NodeInfo* bc_node  = find_input(g, node.id, 0);
    const NodeInfo* rgh_node = find_input(g, node.id, 1);
    const NodeInfo* met_node = find_input(g, node.id, 2);

    // Build uniform declarations and JSON values
    std::ostringstream udecl, ujson, frag;

    // --- baseColor uniform ---
    std::string bc_uname = pfx + "baseColor";
    if (bc_node && bc_node->type == "rgb") {
        baseColor = json_vec3_val(bc_node->props_json, "color", baseColor);
    }
    udecl << "uniform vec3  " << bc_uname << ";\n";
    ujson << "\"" << bc_uname << "\":[" << baseColor[0] << "," << baseColor[1] << "," << baseColor[2] << "]";

    // --- roughness uniform ---
    std::string rgh_uname = pfx + "roughness";
    if (rgh_node && rgh_node->type == "value") {
        roughness = jnum(rgh_node->props_json, "value", roughness);
    }
    udecl << "uniform float " << rgh_uname << ";\n";
    ujson << ",\"" << rgh_uname << "\":" << roughness;

    // --- metalness uniform ---
    std::string met_uname = pfx + "metalness";
    if (met_node && met_node->type == "value") {
        metalness = jnum(met_node->props_json, "value", metalness);
    }
    udecl << "uniform float " << met_uname << ";\n";
    ujson << ",\"" << met_uname << "\":" << metalness;

    // --- alpha uniform ---
    std::string alp_uname = pfx + "alpha";
    udecl << "uniform float " << alp_uname << ";\n";
    ujson << ",\"" << alp_uname << "\":" << alpha_val;

    // --- emission color uniform ---
    std::string emc_uname = pfx + "emissionColor";
    udecl << "uniform vec3  " << emc_uname << ";\n";
    ujson << ",\"" << emc_uname << "\":[" << emColor[0] << "," << emColor[1] << "," << emColor[2] << "]";

    // --- emission strength uniform ---
    std::string ems_uname = pfx + "emissionStrength";
    udecl << "uniform float " << ems_uname << ";\n";
    ujson << ",\"" << ems_uname << "\":" << emStrength;

    // Build fragment body — substitute actual uniform names into PBR template
    // (We inline-expand because the prefix may differ for mix sub-shaders)
    frag << R"(
{
    vec3 N = normalize(vNormal);
    vec3 V = normalize(vViewDir);
    const float PI = 3.14159265359;

    vec3 L0  = normalize(vec3(2.0,4.0,3.0));
    vec3 L1  = normalize(vec3(-3.0,1.0,-2.0));
    vec3 Lc0 = vec3(3.0,2.8,2.5);
    vec3 Lc1 = vec3(0.4,0.45,0.6);

    vec3 albedo = )" << bc_uname << R"(;
    float rough = clamp()" << rgh_uname << R"(, 0.04, 1.0);
    float metal = clamp()" << met_uname << R"(, 0.0,  1.0);
    vec3  F0    = mix(vec3(0.04), albedo, metal);

    // GGX micro-facet evaluation for a single light
    #define PBR_SHADE(Ldir, Lcol) {                                         \
        vec3  H     = normalize(V + (Ldir));                                \
        float NdotL = max(dot(N,(Ldir)),0.0);                               \
        float NdotV = max(dot(N,V),1e-4);                                   \
        float NdotH = max(dot(N,H),0.0);                                    \
        float VdotH = max(dot(V,H),0.0);                                    \
        float a2    = rough*rough; a2=a2*a2;                                \
        float denom = (NdotH*NdotH*(a2-1.0)+1.0);                          \
        float D     = a2/(PI*denom*denom);                                  \
        float kr    = (rough+1.0); kr=kr*kr/8.0;                           \
        float Gv    = NdotV/(NdotV*(1.0-kr)+kr);                           \
        float Gl    = NdotL/(NdotL*(1.0-kr)+kr);                           \
        vec3  F     = F0+(1.0-F0)*pow(clamp(1.0-VdotH,0.0,1.0),5.0);     \
        vec3  kd    = (1.0-F)*(1.0-metal);                                 \
        vec3  spec  = (D*Gv*Gl*F)/max(4.0*NdotV*NdotL,1e-4);             \
        Lo += (kd*albedo/PI + spec)*(Lcol)*NdotL;                          \
    }

    vec3 Lo = vec3(0.0);
    PBR_SHADE(L0, Lc0)
    PBR_SHADE(L1, Lc1)
    #undef PBR_SHADE

    vec3 F_amb = F0+(1.0-F0)*pow(clamp(1.0-max(dot(N,V),0.0),0.0,1.0),5.0);
    vec3 color = Lo
               + vec3(0.03)*albedo
               + F_amb*(1.0-rough)*0.15*albedo
               + )" << emc_uname << R"( * )" << ems_uname << R"(;
    // ACES tone map
    color = clamp((color*(2.51*color+0.03))/(color*(2.43*color+0.59)+0.14),0.0,1.0);
    // Gamma
    color = pow(color, vec3(1.0/2.2));
    gl_FragColor = vec4(color, )" << alp_uname << R"();
}
)";

    ShaderResult sr;
    sr.frag         = frag.str();
    sr.uniforms_decl = udecl.str();
    sr.uniforms_json = ujson.str();
    return sr;
}

// ---------- Emission node ---------------------------------------------------

static ShaderResult compile_emission(const Graph& /*g*/, const NodeInfo& node,
                                      const std::string& pfx) {
    const std::string& props = node.props_json;
    std::array<double,3> col = json_vec3_val(props, "color",    {1.0, 1.0, 1.0});
    double strength          = jnum(props, "strength", 1.0);

    std::string cn = pfx + "emColor";
    std::string sn = pfx + "emStrength";

    std::ostringstream udecl, ujson, frag;
    udecl << "uniform vec3  " << cn << ";\n"
          << "uniform float " << sn << ";\n";
    ujson << "\"" << cn << "\":[" << col[0] << "," << col[1] << "," << col[2] << "]"
          << ",\"" << sn << "\":" << strength;

    frag << "{\n"
         << "    vec3 c = " << cn << " * " << sn << ";\n"
         << "    // Gamma\n"
         << "    c = pow(clamp(c,0.0,1.0), vec3(1.0/2.2));\n"
         << "    gl_FragColor = vec4(c, 1.0);\n"
         << "}\n";

    ShaderResult sr;
    sr.frag          = frag.str();
    sr.uniforms_decl = udecl.str();
    sr.uniforms_json = ujson.str();
    return sr;
}

// ---------- Mix node --------------------------------------------------------

static ShaderResult compile_mix(const Graph& g, const NodeInfo& node,
                                 const std::string& pfx) {
    const std::string& props = node.props_json;
    double factor = jnum(props, "factor", 0.5);

    std::string fn = pfx + "mixFactor";

    // Inputs: port 0 = shader1, port 1 = shader2
    const NodeInfo* s1 = find_input(g, node.id, 0);
    const NodeInfo* s2 = find_input(g, node.id, 1);

    // Compile sub-shaders with unique prefixes
    ShaderResult r1, r2;
    bool have_s1 = false, have_s2 = false;

    if (s1) { r1 = compile_node(g, *s1, pfx + "s1_"); have_s1 = true; }
    if (s2) { r2 = compile_node(g, *s2, pfx + "s2_"); have_s2 = true; }

    std::ostringstream udecl, ujson, frag;

    udecl << "uniform float " << fn << ";\n";
    ujson << "\"" << fn << "\":" << factor;

    if (have_s1) {
        udecl << r1.uniforms_decl;
        if (!r1.uniforms_json.empty()) ujson << "," << r1.uniforms_json;
    }
    if (have_s2) {
        udecl << r2.uniforms_decl;
        if (!r2.uniforms_json.empty()) ujson << "," << r2.uniforms_json;
    }

    // We evaluate both sub-shaders into temporary gl_FragColor values and mix
    frag << "{\n"
         << "    vec4 mix_c1 = vec4(0.5,0.0,0.5,1.0);\n"
         << "    vec4 mix_c2 = vec4(0.5,0.0,0.5,1.0);\n";
    if (have_s1) {
        frag << "    // Shader 1\n"
             << "    gl_FragColor = vec4(0.0);\n"
             << "    " << r1.frag << "\n"
             << "    mix_c1 = gl_FragColor;\n";
    }
    if (have_s2) {
        frag << "    // Shader 2\n"
             << "    gl_FragColor = vec4(0.0);\n"
             << "    " << r2.frag << "\n"
             << "    mix_c2 = gl_FragColor;\n";
    }
    frag << "    gl_FragColor = mix(mix_c1, mix_c2, " << fn << ");\n"
         << "}\n";

    ShaderResult sr;
    sr.frag          = frag.str();
    sr.uniforms_decl = udecl.str();
    sr.uniforms_json = ujson.str();
    return sr;
}

// ---------- Default (purple Phong) ------------------------------------------

static ShaderResult compile_default() {
    ShaderResult sr;
    sr.uniforms_decl = "";
    sr.uniforms_json = "\"defaultColor\":[0.5,0.0,0.5],\"defaultShininess\":32.0";
    sr.frag = R"({
    vec3 N = normalize(vNormal);
    vec3 V = normalize(vViewDir);
    vec3 L = normalize(vec3(2.0,4.0,3.0));
    vec3 albedo = vec3(0.5, 0.0, 0.5);
    float diff = max(dot(N,L),0.0);
    vec3 R = reflect(-L, N);
    float spec = pow(max(dot(V,R),0.0), 32.0);
    vec3 color = albedo*(0.1+diff*0.9) + vec3(0.5)*spec;
    color = pow(clamp(color,0.0,1.0), vec3(1.0/2.2));
    gl_FragColor = vec4(color, 1.0);
}
)";
    return sr;
}

// ---------- Dispatch --------------------------------------------------------

static ShaderResult compile_node(const Graph& g, const NodeInfo& node,
                                  const std::string& prefix) {
    if (node.type == "principled") return compile_principled(g, node, prefix);
    if (node.type == "emission")   return compile_emission(g, node, prefix);
    if (node.type == "mix")        return compile_mix(g, node, prefix);
    return compile_default();
}

// =============================================================================
// Public compile() entry point
// =============================================================================

std::string compile(const std::string& graph_json) {
    // --- Parse ---
    Graph g;
    try {
        g = parse_graph(graph_json);
    } catch (...) {
        return "{\"vertex\":\"\",\"fragment\":\"\",\"uniforms\":{},\"error\":\"Failed to parse graph JSON\"}";
    }

    // --- Find output node ---
    const NodeInfo* output_node = nullptr;
    for (const auto& kv : g.nodes) {
        if (kv.second.type == "output") { output_node = &kv.second; break; }
    }
    if (!output_node) {
        return "{\"vertex\":\"\",\"fragment\":\"\",\"uniforms\":{},\"error\":\"No output node found\"}";
    }

    // --- Find node connected to output port 0 (Surface) ---
    const NodeInfo* surface_node = find_input(g, output_node->id, 0);

    ShaderResult sr;
    if (surface_node) {
        sr = compile_node(g, *surface_node, "");
    } else {
        sr = compile_default();
    }

    // --- Build full fragment shader source ---
    std::ostringstream frag_full;
    if (!sr.uniforms_decl.empty())
        frag_full << sr.uniforms_decl << "\n";
    frag_full << "varying vec3 vNormal;\n"
              << "varying vec3 vWorldPos;\n"
              << "varying vec3 vViewDir;\n"
              << "\n"
              << "void main() " << sr.frag << "\n";

    // --- JSON-escape GLSL strings ---
    auto escape = [](const std::string& s) {
        std::string out;
        out.reserve(s.size() + 32);
        for (unsigned char c : s) {
            if      (c == '"')  out += "\\\"";
            else if (c == '\\') out += "\\\\";
            else if (c == '\n') out += "\\n";
            else if (c == '\r') out += "\\r";
            else if (c == '\t') out += "\\t";
            else if (c < 0x20) {
                out += "\\u00";
                out += "0123456789abcdef"[c >> 4];
                out += "0123456789abcdef"[c & 0xf];
            } else {
                out += static_cast<char>(c);
            }
        }
        return out;
    };

    std::string vert_escaped = escape(VERTEX_GLSL);
    std::string frag_escaped = escape(frag_full.str());

    std::ostringstream result;
    result << "{"
           << "\"vertex\":\""   << vert_escaped << "\","
           << "\"fragment\":\"" << frag_escaped << "\","
           << "\"uniforms\":{"  << sr.uniforms_json << "},"
           << "\"error\":\"\""
           << "}";

    return result.str();
}

} // namespace node_compiler
