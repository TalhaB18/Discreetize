#include "shader_library.hpp"
#include <sstream>

namespace shaders {

// ---------------------------------------------------------------------------
// PBR — Physically-based shading
// ---------------------------------------------------------------------------

static const char PBR_VERT[] = R"glsl(
varying vec3 vNormal;
varying vec3 vViewPos;
varying vec3 vWorldPos;

void main() {
    vec4 mvPos = modelViewMatrix * vec4(position, 1.0);
    vViewPos  = -mvPos.xyz;
    vWorldPos = (modelMatrix * vec4(position, 1.0)).xyz;
    vNormal   = normalize(normalMatrix * normal);
    gl_Position = projectionMatrix * mvPos;
}
)glsl";

static const char PBR_FRAG[] = R"glsl(
varying vec3 vNormal;
varying vec3 vViewPos;
varying vec3 vWorldPos;

uniform vec3  uColor;
uniform float uRoughness;

const vec3  LIGHT_DIR  = normalize(vec3(1.0, 2.0, 3.0));
const vec3  LIGHT_COL  = vec3(1.0, 0.98, 0.92);
const float METALNESS  = 0.1;

float DistributionGGX(vec3 N, vec3 H, float rough) {
    float a    = rough * rough;
    float a2   = a * a;
    float NdotH = max(dot(N, H), 0.0);
    float denom = (NdotH * NdotH * (a2 - 1.0) + 1.0);
    return a2 / (3.14159265 * denom * denom + 0.0001);
}

float GeometrySchlick(float NdotV, float rough) {
    float r = rough + 1.0;
    float k = (r * r) / 8.0;
    return NdotV / (NdotV * (1.0 - k) + k);
}

vec3 fresnelSchlick(float cosTheta, vec3 F0) {
    return F0 + (1.0 - F0) * pow(clamp(1.0 - cosTheta, 0.0, 1.0), 5.0);
}

void main() {
    vec3 N = normalize(vNormal);
    vec3 V = normalize(vViewPos);
    vec3 L = LIGHT_DIR;
    vec3 H = normalize(V + L);

    vec3 F0 = mix(vec3(0.04), uColor, METALNESS);

    float NDF = DistributionGGX(N, H, uRoughness);
    float G   = GeometrySchlick(max(dot(N, V), 0.0), uRoughness)
              * GeometrySchlick(max(dot(N, L), 0.0), uRoughness);
    vec3  F   = fresnelSchlick(max(dot(H, V), 0.0), F0);

    float denom    = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0) + 0.0001;
    vec3  specular = (NDF * G * F) / denom;

    vec3 kD = (1.0 - F) * (1.0 - METALNESS);
    vec3 diffuse = kD * uColor / 3.14159265;

    float NdotL = max(dot(N, L), 0.0);
    vec3  Lo    = (diffuse + specular) * LIGHT_COL * NdotL;

    // IBL ambient approximation
    vec3 ambient = vec3(0.03) * uColor + 0.06 * uColor * max(dot(N, vec3(0.0, 1.0, 0.0)), 0.0);

    vec3 color = ambient + Lo;
    // Tone map (Reinhard)
    color = color / (color + vec3(1.0));
    color = pow(color, vec3(1.0 / 2.2));

    gl_FragColor = vec4(color, 1.0);
}
)glsl";

// ---------------------------------------------------------------------------
// Matcap — sphere capture
// ---------------------------------------------------------------------------

static const char MATCAP_VERT[] = R"glsl(
varying vec2 vMatcapUV;

void main() {
    vec3 nView = normalize((modelViewMatrix * vec4(normal, 0.0)).xyz);
    vMatcapUV = nView.xy * 0.5 + 0.5;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
)glsl";

static const char MATCAP_FRAG[] = R"glsl(
varying vec2 vMatcapUV;

void main() {
    // Procedural matcap: gradient from deep purple at bottom to bright white at top
    vec3 purple = vec3(0.42, 0.34, 0.80);
    vec3 midCol = vec3(0.70, 0.60, 0.98);
    vec3 white  = vec3(1.00, 1.00, 1.00);

    float t = vMatcapUV.y;
    vec3 col;
    if (t < 0.5) {
        col = mix(purple, midCol, t * 2.0);
    } else {
        col = mix(midCol, white, (t - 0.5) * 2.0);
    }

    // Specular highlight near top-right
    float spec = pow(max(0.0, vMatcapUV.x * 0.6 + vMatcapUV.y * 0.6 - 0.5), 4.0);
    col = mix(col, vec3(1.0), spec * 0.8);

    // Rim darkening at edges
    float rim = length(vMatcapUV * 2.0 - 1.0);
    col = mix(col, col * 0.3, smoothstep(0.7, 1.0, rim));

    gl_FragColor = vec4(col, 1.0);
}
)glsl";

// ---------------------------------------------------------------------------
// Normal — normal visualization
// ---------------------------------------------------------------------------

static const char NORMAL_VERT[] = R"glsl(
varying vec3 vNormal;

void main() {
    vNormal = normalize((modelMatrix * vec4(normal, 0.0)).xyz);
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
)glsl";

static const char NORMAL_FRAG[] = R"glsl(
varying vec3 vNormal;

void main() {
    vec3 n = normalize(vNormal);
    gl_FragColor = vec4(n * 0.5 + 0.5, 1.0);
}
)glsl";

// ---------------------------------------------------------------------------
// Wireframe — barycentric wireframe overlay
// ---------------------------------------------------------------------------

static const char WIREFRAME_VERT[] = R"glsl(
attribute vec3 barycentric;
varying vec3 vBarycentric;
varying vec3 vNormal;

void main() {
    vBarycentric = barycentric;
    vNormal = normalize(normalMatrix * normal);
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
)glsl";

static const char WIREFRAME_FRAG[] = R"glsl(
varying vec3 vBarycentric;
varying vec3 vNormal;

float edgeFactor() {
    vec3 d = fwidth(vBarycentric);
    vec3 a = smoothstep(vec3(0.0), d * 1.2, vBarycentric);
    return min(min(a.x, a.y), a.z);
}

void main() {
    float e = edgeFactor();
    // White edges, transparent faces
    vec3  edgeCol = vec3(1.0, 1.0, 1.0);
    float alpha   = 1.0 - e;         // edge: alpha=1, face interior: alpha=0
    if (alpha < 0.01) discard;
    gl_FragColor = vec4(edgeCol, alpha);
}
)glsl";

// ---------------------------------------------------------------------------
// X-ray — Fresnel transparency
// ---------------------------------------------------------------------------

static const char XRAY_VERT[] = R"glsl(
varying vec3 vNormal;
varying vec3 vViewDir;

void main() {
    vec4 mvPos = modelViewMatrix * vec4(position, 1.0);
    vNormal  = normalize(normalMatrix * normal);
    vViewDir = normalize(-mvPos.xyz);
    gl_Position = projectionMatrix * mvPos;
}
)glsl";

static const char XRAY_FRAG[] = R"glsl(
varying vec3 vNormal;
varying vec3 vViewDir;

uniform vec3 uColor;

void main() {
    vec3  N   = normalize(vNormal);
    vec3  V   = normalize(vViewDir);
    float rim = pow(1.0 - abs(dot(V, N)), 2.0);
    gl_FragColor = vec4(uColor, rim);
}
)glsl";

// ---------------------------------------------------------------------------
// Depth — linear depth visualization
// ---------------------------------------------------------------------------

static const char DEPTH_VERT[] = R"glsl(
varying float vDepth;

void main() {
    vec4 clip = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    vDepth = clip.z / clip.w;
    gl_Position = clip;
}
)glsl";

static const char DEPTH_FRAG[] = R"glsl(
varying float vDepth;

void main() {
    // Linearize NDC depth [-1,1] -> [0,1], then remap for visibility
    float d = vDepth * 0.5 + 0.5;
    // Apply power curve to spread out near-geometry contrast
    d = pow(d, 0.4);
    gl_FragColor = vec4(vec3(d), 1.0);
}
)glsl";

// ---------------------------------------------------------------------------
// Toon — cel shading with hard outline
// ---------------------------------------------------------------------------

static const char TOON_VERT[] = R"glsl(
varying vec3 vNormal;
varying vec3 vViewDir;

void main() {
    vec4 mvPos = modelViewMatrix * vec4(position, 1.0);
    vNormal  = normalize(normalMatrix * normal);
    vViewDir = normalize(-mvPos.xyz);
    gl_Position = projectionMatrix * mvPos;
}
)glsl";

static const char TOON_FRAG[] = R"glsl(
varying vec3 vNormal;
varying vec3 vViewDir;

uniform vec3 uColor;

void main() {
    vec3  N     = normalize(vNormal);
    vec3  V     = normalize(vViewDir);
    vec3  L     = normalize(vec3(1.0, 2.0, 3.0));
    float NdotL = max(dot(N, L), 0.0);

    // Quantize into 4 bands
    float bands  = 4.0;
    float band   = floor(NdotL * bands) / bands;

    // Hard outline: near-silhouette pixels become black
    float edge   = dot(N, V);
    if (edge < 0.25) {
        gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);
        return;
    }

    vec3 shadedCol = uColor * (0.2 + 0.8 * band);
    gl_FragColor   = vec4(shadedCol, 1.0);
}
)glsl";

// ---------------------------------------------------------------------------
// Registry
// ---------------------------------------------------------------------------

static std::map<std::string, ShaderDef> build_registry() {
    std::map<std::string, ShaderDef> m;

    m["pbr"] = {
        "pbr",
        "Physically-based shading with GGX specular and IBL-approximate ambient",
        PBR_VERT, PBR_FRAG,
        "surface"
    };
    m["matcap"] = {
        "matcap",
        "Matcap sphere capture with procedural purple-to-white gradient",
        MATCAP_VERT, MATCAP_FRAG,
        "surface"
    };
    m["normal"] = {
        "normal",
        "World-space normal visualization (RGB = normal.xyz * 0.5 + 0.5)",
        NORMAL_VERT, NORMAL_FRAG,
        "debug"
    };
    m["wireframe"] = {
        "wireframe",
        "Barycentric wireframe overlay — white edges on transparent faces",
        WIREFRAME_VERT, WIREFRAME_FRAG,
        "debug"
    };
    m["xray"] = {
        "xray",
        "Fresnel X-ray transparency — silhouette edges opaque, interior transparent",
        XRAY_VERT, XRAY_FRAG,
        "debug"
    };
    m["depth"] = {
        "depth",
        "Linearized depth visualization as greyscale gradient",
        DEPTH_VERT, DEPTH_FRAG,
        "debug"
    };
    m["toon"] = {
        "toon",
        "Cel shading with 4-band diffuse quantization and hard silhouette outline",
        TOON_VERT, TOON_FRAG,
        "stylized"
    };

    return m;
}

const std::map<std::string, ShaderDef>& get_all() {
    static const auto registry = build_registry();
    return registry;
}

const ShaderDef* get(const std::string& name) {
    const auto& m = get_all();
    auto it = m.find(name);
    return (it != m.end()) ? &it->second : nullptr;
}

// JSON-escape a GLSL source string: escape backslashes, quotes, and newlines.
static std::string json_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 64);
    for (char c : s) {
        if      (c == '\\') out += "\\\\";
        else if (c == '"')  out += "\\\"";
        else if (c == '\n') out += "\\n";
        else if (c == '\r') out += "\\r";
        else if (c == '\t') out += "\\t";
        else                out += c;
    }
    return out;
}

std::string list_json() {
    // Emit in a stable order matching definition order above.
    static const std::vector<std::string> order = {
        "pbr", "matcap", "normal", "wireframe", "xray", "depth", "toon"
    };
    const auto& m = get_all();
    std::ostringstream o;
    o << "[";
    bool first = true;
    for (const auto& key : order) {
        auto it = m.find(key);
        if (it == m.end()) continue;
        if (!first) o << ",";
        first = false;
        const auto& s = it->second;
        o << "{\"name\":\"" << s.name << "\","
          << "\"description\":\"" << json_escape(s.description) << "\","
          << "\"category\":\"" << s.category << "\"}";
    }
    o << "]";
    return o.str();
}

std::string shader_json(const std::string& name) {
    const ShaderDef* s = get(name);
    if (!s) return "";
    std::ostringstream o;
    o << "{"
      << "\"name\":\""        << s->name                      << "\","
      << "\"description\":\"" << json_escape(s->description)  << "\","
      << "\"category\":\""    << s->category                  << "\","
      << "\"vert\":\""        << json_escape(s->vert_src)     << "\","
      << "\"frag\":\""        << json_escape(s->frag_src)     << "\""
      << "}";
    return o.str();
}

} // namespace shaders
