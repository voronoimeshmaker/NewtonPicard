#pragma once

#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiDiagram2D.hpp>

#include "gsmaker/core/Config.hpp"
#include "gsmaker/core/MeshAdapter.hpp"

namespace gss {

/// @brief Constructs solver-facing meshes from mesh-only configuration.
///
/// MeshFactory owns mesh-generation policy and the interaction with
/// VoronoiMeshMaker. MeshAdapter remains responsible for storing and
/// validating the solver-facing mesh representation.
///
/// The explicit diagram-producing path exists so numerical modules that depend
/// on the same physical Voronoi geometry (notably VoroGradToolkit) can share
/// one generated diagram instead of regenerating an equivalent mesh.
class MeshFactory final {
public:
    using VoronoiDiagram2D =
        vmm::vd2d::ClippedVoronoiDiagram2D;

    [[nodiscard]]
    static MeshAdapter create(
        const MeshConfig& config);

    [[nodiscard]]
    static MeshAdapter create(
        const VoronoiDiagram2D& diagram);

    [[nodiscard]]
    static VoronoiDiagram2D create_voronoi_diagram(
        const MeshConfig& config);

private:
    [[nodiscard]]
    static MeshAdapter create_uniform_internal(
        const MeshConfig& config);
};

} // namespace gss


