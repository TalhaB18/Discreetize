#pragma once
#include <string>

namespace node_compiler {

// Compile a node graph JSON into GLSL shaders.
// Input JSON format:
// {
//   "nodes": [
//     {"id":0,"type":"output"},
//     {"id":1,"type":"principled","props":{"baseColor":[0.8,0.2,0.2],"roughness":0.5,"metalness":0.0,"alpha":1.0,"emissionColor":[0,0,0],"emissionStrength":0.0}},
//     {"id":2,"type":"rgb","props":{"color":[0.8,0.2,0.2]}},
//     {"id":3,"type":"value","props":{"value":0.5}},
//     {"id":4,"type":"emission","props":{"color":[1,0.5,0],"strength":2.0}},
//     {"id":5,"type":"mix","props":{"factor":0.5}}
//   ],
//   "connections": [
//     {"from":1,"fromPort":0,"to":0,"toPort":0},
//     {"from":2,"fromPort":0,"to":1,"toPort":0}
//   ]
// }
//
// Returns JSON: {"vertex":"...glsl...","fragment":"...glsl...","uniforms":{...},"error":""}
std::string compile(const std::string& graph_json);

} // namespace node_compiler
