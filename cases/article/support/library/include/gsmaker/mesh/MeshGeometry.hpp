#pragma once

#include <span>
#include <vector>

#include "gsmaker/core/Types.hpp"
#include "gsmaker/mesh/BoundaryFaceGeometry.hpp"

namespace gss {

/// @brief Pure geometric data attached to the finite-volume mesh.
class MeshGeometry {
public:
    MeshGeometry() = default;
    explicit MeshGeometry(Index n_cells);
    explicit MeshGeometry(std::vector<Real> cell_volume);

    void resize(Index n_cells);
    void set_cell_volume(std::vector<Real> cell_volume);
    void set_boundary_faces(std::vector<BoundaryFaceGeometry> boundary_faces);

    [[nodiscard]] Index n_cells() const noexcept { return n_cells_; }

    [[nodiscard]] std::span<Real> cell_volume() noexcept { return cell_volume_; }
    [[nodiscard]] std::span<const Real> cell_volume() const noexcept { return cell_volume_; }

    [[nodiscard]] std::span<BoundaryFaceGeometry> boundary_faces() noexcept {
        return boundary_faces_;
    }
    [[nodiscard]] std::span<const BoundaryFaceGeometry> boundary_faces() const noexcept {
        return boundary_faces_;
    }

    [[nodiscard]] Real total_volume() const noexcept;
    void validate() const;

private:
    Index n_cells_{0};
    std::vector<Real> cell_volume_{};
    std::vector<BoundaryFaceGeometry> boundary_faces_{};
};

} // namespace gss


