#include "gsmaker/mesh/MeshGeometry.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace gss {
namespace {

[[nodiscard]] std::size_t checked_size(Index n_cells) {
    if (n_cells < Index{0}) {
        throw std::invalid_argument("MeshGeometry cell count must be non-negative");
    }
    return static_cast<std::size_t>(n_cells);
}

[[nodiscard]] bool finite(Real value) noexcept {
    return std::isfinite(static_cast<double>(value));
}

} // namespace

MeshGeometry::MeshGeometry(Index n_cells) { resize(n_cells); }

MeshGeometry::MeshGeometry(std::vector<Real> cell_volume) {
    set_cell_volume(std::move(cell_volume));
}

void MeshGeometry::resize(Index n_cells) {
    const auto size = checked_size(n_cells);
    n_cells_ = n_cells;
    cell_volume_.assign(size, Real{0});
    boundary_faces_.clear();
}

void MeshGeometry::set_cell_volume(std::vector<Real> cell_volume) {
    n_cells_ = static_cast<Index>(cell_volume.size());
    cell_volume_ = std::move(cell_volume);
}

void MeshGeometry::set_boundary_faces(
    std::vector<BoundaryFaceGeometry> boundary_faces) {
    boundary_faces_ = std::move(boundary_faces);
}

Real MeshGeometry::total_volume() const noexcept {
    return std::accumulate(cell_volume_.begin(), cell_volume_.end(), Real{0});
}

void MeshGeometry::validate() const {
    const auto expected = checked_size(n_cells_);
    if (cell_volume_.size() != expected) {
        throw std::runtime_error("MeshGeometry cell-volume array does not match cell count");
    }

    for (const Real volume : cell_volume_) {
        if (!finite(volume) || !(volume > Real{0})) {
            throw std::runtime_error("MeshGeometry cell volumes must be finite and strictly positive");
        }
    }

    constexpr double normal_tolerance = 1.0e-6;
    for (const auto& face : boundary_faces_) {
        if (!finite(face.x) || !finite(face.y) ||
            !finite(face.normal_x) || !finite(face.normal_y) ||
            !finite(face.measure) || !finite(face.perpendicular_distance)) {
            throw std::runtime_error("Boundary-face geometry contains a non-finite value");
        }
        if (!(face.measure > Real{0}) || !(face.perpendicular_distance > Real{0})) {
            throw std::runtime_error("Boundary-face measure and perpendicular distance must be positive");
        }
        const double normal_norm = std::hypot(
            static_cast<double>(face.normal_x),
            static_cast<double>(face.normal_y));
        if (std::abs(normal_norm - 1.0) > normal_tolerance) {
            throw std::runtime_error("Boundary-face outward normal must have unit length");
        }
    }
}

} // namespace gss
