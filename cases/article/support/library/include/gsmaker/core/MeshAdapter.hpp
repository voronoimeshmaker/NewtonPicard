#pragma once

#include <vector>

#include "gsmaker/core/Config.hpp"
#include "gsmaker/core/Types.hpp"
#include "gsmaker/mesh/BoundaryFaceGeometry.hpp"
#include "gsmaker/mesh/MeshGeometry.hpp"


namespace gss {

/// @brief Solver-facing representation of mesh topology and pure geometry.
class MeshAdapter {
public:
    MeshAdapter() = default;

    MeshAdapter(
        MeshTopology topology,
        std::vector<Real> volumes,
        std::vector<BoundaryFaceGeometry> boundary_faces = {});

    [[nodiscard]] static MeshAdapter uniform_from_config(const Config& config);
    [[nodiscard]] static MeshAdapter from_config(const Config& config);

    [[nodiscard]] Index n_cells() const noexcept { return geometry_.n_cells(); }
    [[nodiscard]] const MeshTopology& topology() const noexcept { return topology_; }

    [[nodiscard]] MeshGeometry& geometry() noexcept { return geometry_; }
    [[nodiscard]] const MeshGeometry& geometry() const noexcept { return geometry_; }

    [[nodiscard]] CellData& data() noexcept { return data_; }
    [[nodiscard]] const CellData& data() const noexcept { return data_; }

    [[nodiscard]] Real h_eff() const noexcept;
    [[nodiscard]] Real total_volume() const noexcept;
    void validate() const;

private:
    MeshTopology topology_{};
    MeshGeometry geometry_{};
    CellData data_{};
};

} // namespace gss
