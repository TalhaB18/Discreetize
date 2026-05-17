#pragma once
#include <string>
#include <vector>
#include <map>

namespace shaders {

struct ShaderDef {
    std::string name;         // "pbr", "matcap", "normal", "wireframe", "xray", "depth", "toon"
    std::string description;
    std::string vert_src;     // GLSL vertex shader source
    std::string frag_src;     // GLSL fragment shader source
    std::string category;     // "surface", "debug", "stylized"
};

const std::map<std::string, ShaderDef>& get_all();
const ShaderDef* get(const std::string& name);
std::string list_json();                          // JSON array of {name, description, category}
std::string shader_json(const std::string& name); // JSON {name, vert, frag, ...}

} // namespace shaders
