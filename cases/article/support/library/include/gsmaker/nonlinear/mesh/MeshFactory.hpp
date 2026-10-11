#pragma once

#include "gsmaker/core/Config.hpp"
#include "gsmaker/core/MeshAdapter.hpp"

namespace gss {

/// @brief Constructs solver-facing meshes from mesh-only configuration.
///
/// MeshFactory owns mesh-generation policy and the interaction with
/// VoronoiMeshMaker. MeshAdapter remains responsible for storing and
/// validating the solver-facing mesh representation.
class MeshFactory final {
public:
    [[nodiscard]]
    static MeshAdapter create(
        const MeshConfig& config);

private:
    [[nodiscard]]
    static MeshAdapter create_uniform_internal(
        const MeshConfig& config);
};

} // namespace gss
