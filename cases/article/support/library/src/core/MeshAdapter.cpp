#include "gsmaker/core/MeshAdapter.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "gsmaker/mesh/MeshFactory.hpp"

namespace gss {

MeshAdapter::MeshAdapter(
    MeshTopology topology,
    std::vector<Real> volumes,
    std::vector<BoundaryFaceGeometry> boundary_faces)
    : topology_(std::move(topology)) {

    geometry_.set_cell_volume(std::move(volumes));
    geometry_.set_boundary_faces(std::move(boundary_faces));
    data_.resize(geometry_.n_cells());
    validate();
}

MeshAdapter MeshAdapter::uniform_from_config(const Config& config) {
    MeshConfig mesh_config = config.mesh;
    mesh_config.type = "uniform_internal";
    return MeshFactory::create(mesh_config);
}

MeshAdapter MeshAdapter::from_config(const Config& config) {
    return MeshFactory::create(config.mesh);
}

Real MeshAdapter::h_eff() const noexcept {
    if (n_cells() <= 0) return Real{0};
    return std::sqrt(total_volume() / static_cast<Real>(n_cells()));
}

Real MeshAdapter::total_volume() const noexcept {
    return geometry_.total_volume();
}

void MeshAdapter::validate() const {
    geometry_.validate();

    const auto n = static_cast<std::size_t>(geometry_.n_cells());
    if (topology_.neighbor_start.size() != n + 1) {
        throw std::runtime_error("mesh CSR neighbor_start must have n_cells + 1 entries");
    }
    if (topology_.coord_x.size() != n || topology_.coord_y.size() != n) {
        throw std::runtime_error("mesh coordinate arrays must match n_cells");
    }
    if (topology_.is_boundary_cell.size() != n) {
        throw std::runtime_error("mesh boundary-cell array must match n_cells");
    }
    if (topology_.neighbor_list.size() != topology_.interface_area.size() ||
        topology_.neighbor_list.size() != topology_.distance.size()) {
        throw std::runtime_error("mesh CSR edge arrays must have equal length");
    }
    if (topology_.neighbor_start.empty() ||
        topology_.neighbor_start.front() != Index{0} ||
        topology_.neighbor_start.back() != static_cast<Index>(topology_.neighbor_list.size())) {
        throw std::runtime_error("mesh CSR offsets are inconsistent");
    }

    for (std::size_t cell = 0; cell < n; ++cell) {
        const Index begin = topology_.neighbor_start[cell];
        const Index end = topology_.neighbor_start[cell + 1];
        if (begin < 0 || end < begin ||
            end > static_cast<Index>(topology_.neighbor_list.size())) {
            throw std::runtime_error("mesh CSR offsets are invalid");
        }
    }

    for (std::size_t edge = 0; edge < topology_.neighbor_list.size(); ++edge) {
        const Index neighbour = topology_.neighbor_list[edge];
        if (neighbour < 0 || neighbour >= static_cast<Index>(n)) {
            throw std::runtime_error("mesh contains an out-of-range neighbour");
        }
        if (!(topology_.interface_area[edge] > Real{0}) ||
            !(topology_.distance[edge] > Real{0}) ||
            !std::isfinite(static_cast<double>(topology_.interface_area[edge])) ||
            !std::isfinite(static_cast<double>(topology_.distance[edge]))) {
            throw std::runtime_error("interface area and distance must be finite and positive");
        }
    }

    for (const auto& face : geometry_.boundary_faces()) {
        if (face.owner_cell < 0 || face.owner_cell >= static_cast<Index>(n)) {
            throw std::runtime_error("boundary face has an out-of-range owner cell");
        }
        if (topology_.is_boundary_cell[static_cast<std::size_t>(face.owner_cell)] == 0U) {
            throw std::runtime_error("boundary face owner is not marked as a boundary cell");
        }
    }
}

} // namespace gss



