#pragma once

#include <stdexcept>

#include "gss/core/Config.hpp"
#include "gss/core/Types.hpp"

namespace gss {

/// @brief Lightweight adapter exposing VoronoiMeshMaker data as solver spans.
class MeshAdapter {
public:
    /// @brief Construct an empty adapter.
    MeshAdapter() = default;

    /// @brief Construct from already exported Voronoi topology and volumes.
    MeshAdapter(MeshTopology topology, std::vector<Real> volumes);

    /// @brief Build a deterministic uniform Voronoi-like Cartesian control-volume mesh for tests and smoke runs.
    [[nodiscard]] static MeshAdapter uniform_from_config(const Config& config);

    /// @brief Build a clipped Voronoi mesh through VoronoiMeshMaker from Config.
    [[nodiscard]] static MeshAdapter from_config(const Config& config);

    /// @brief Number of control volumes.
    [[nodiscard]] Index n_cells() const noexcept { return static_cast<Index>(data_.n_cells); }

    /// @brief Read-only topology.
    [[nodiscard]] const MeshTopology& topology() const noexcept { return topology_; }

    /// @brief Cell renumbering permutation (new slot -> original physical cell).
    /// Empty when the mesh uses the natural ordering.
    [[nodiscard]] const std::vector<Index>& ordering() const noexcept { return ordering_; }

    /// @brief Record the renumbering permutation used to build this mesh.
    void set_ordering(std::vector<Index> perm) { ordering_ = std::move(perm); }

    /// @brief Mutable hot cell data.
    [[nodiscard]] CellData& data() noexcept { return data_; }

    /// @brief Read-only hot cell data.
    [[nodiscard]] const CellData& data() const noexcept { return data_; }

    /// @brief Effective mesh size sqrt(total volume / number of cells).
    [[nodiscard]] Real h_eff() const noexcept;

    /// @brief Total cell area.
    [[nodiscard]] Real total_volume() const noexcept;

    /// @brief Validate CSR sizes and metric positivity.
    void validate() const;

private:
    MeshTopology topology_{};
    CellData data_{};
    std::vector<Index> ordering_{};
};

} // namespace gss
