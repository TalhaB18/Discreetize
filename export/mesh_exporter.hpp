#pragma once
#include "../geometry/types.hpp"
#include "../meshing/mesh_types.hpp"
#include <memory>
#include <string>
#include <vector>
#include <array>

namespace cfd {

struct CurvatureResult;

class MeshExporter {
public:
    virtual ~MeshExporter() = default;
    virtual void write(const SurfaceManifold& m, const std::string& path) const = 0;
    virtual std::string extension() const = 0;
};

class STLExporter : public MeshExporter {
public:
    enum class Mode { Binary, ASCII };
    explicit STLExporter(Mode mode = Mode::Binary) : mode_(mode) {}
    void write(const SurfaceManifold& m, const std::string& path) const override;
    std::string extension() const override;
private:
    Mode mode_;
    void writeBinary(const SurfaceManifold& m, const std::string& path) const;
    void writeASCII (const SurfaceManifold& m, const std::string& path) const;
};

class OBJExporter : public MeshExporter {
public:
    explicit OBJExporter(const std::vector<std::vector<uint32_t>>* groups = nullptr)
        : groups_(groups) {}
    void write(const SurfaceManifold& m, const std::string& path) const override;
    std::string extension() const override;
private:
    const std::vector<std::vector<uint32_t>>* groups_;
};

class PLYExporter : public MeshExporter {
public:
    explicit PLYExporter(const CurvatureResult* curvature = nullptr)
        : curvature_(curvature) {}
    void write(const SurfaceManifold& m, const std::string& path) const override;
    std::string extension() const override;
private:
    const CurvatureResult* curvature_;
};

class STEPExporter : public MeshExporter {
public:
    void write(const SurfaceManifold& m, const std::string& path) const override;
    std::string extension() const override;
};

// VTU = VTK XML UnstructuredGrid (readable by ParaView)
class VTUExporter : public MeshExporter {
public:
    // Write surface triangles (VTK type 5)
    void write(const SurfaceManifold& m, const std::string& path) const override;
    // Write hybrid volume mesh (tets = type 10, prisms = type 13)
    // Cell data arrays: "cell_type" (0=tet,1=prism) and "layer" for prisms.
    void writeVolume(const VolumeMesh& vol, const std::string& path) const;
    std::string extension() const override { return ".vtu"; }
};

std::unique_ptr<MeshExporter> makeExporter(
    const std::string& format,
    const std::vector<std::vector<uint32_t>>* groups   = nullptr,
    const CurvatureResult*                    curvature = nullptr);

} // namespace cfd
