#include "gss/core/CellPermutation.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace gss {
namespace {

[[nodiscard]] std::vector<Index> make_old_to_new(std::span<const Index> new_to_old) {
    validate_cell_permutation(new_to_old);

    std::vector<Index> old_to_new(new_to_old.size());
    for (Index new_index = 0; new_index < static_cast<Index>(new_to_old.size()); ++new_index) {
        const Index old_index = new_to_old[static_cast<std::size_t>(new_index)];
        old_to_new[static_cast<std::size_t>(old_index)] = new_index;
    }
    return old_to_new;
}

} // namespace

CellPermutation::CellPermutation(std::vector<Index> new_to_old, std::vector<Index> old_to_new)
    : new_to_old_(std::move(new_to_old)), old_to_new_(std::move(old_to_new)) {}

CellPermutation CellPermutation::identity(Index n) {
    if (n < 0) {
        throw std::runtime_error("cell permutation size must be non-negative");
    }
    std::vector<Index> mapping(static_cast<std::size_t>(n));
    std::iota(mapping.begin(), mapping.end(), Index{0});
    return from_new_to_old(std::move(mapping));
}

CellPermutation CellPermutation::from_new_to_old(std::vector<Index> new_to_old) {
    auto old_to_new = make_old_to_new(new_to_old);
    return CellPermutation(std::move(new_to_old), std::move(old_to_new));
}

bool CellPermutation::is_identity() const noexcept {
    for (Index i = 0; i < static_cast<Index>(new_to_old_.size()); ++i) {
        if (new_to_old_[static_cast<std::size_t>(i)] != i) {
            return false;
        }
    }
    return true;
}

CellPermutation CellPermutation::inverse() const {
    return from_new_to_old(old_to_new_);
}

CellPermutation CellPermutation::then(const CellPermutation& next) const {
    if (size() != next.size()) {
        throw std::runtime_error("cannot compose cell permutations with different sizes");
    }

    std::vector<Index> composed(new_to_old_.size());
    for (Index i = 0; i < static_cast<Index>(composed.size()); ++i) {
        const Index intermediate = next.new_to_old_[static_cast<std::size_t>(i)];
        composed[static_cast<std::size_t>(i)] = new_to_old_[static_cast<std::size_t>(intermediate)];
    }
    return from_new_to_old(std::move(composed));
}

void validate_cell_permutation(std::span<const Index> new_to_old) {
    std::vector<unsigned char> seen(new_to_old.size(), 0U);
    for (const Index old_index : new_to_old) {
        if (old_index < 0 || old_index >= static_cast<Index>(new_to_old.size())) {
            throw std::runtime_error("cell permutation contains an out-of-range index");
        }
        auto& slot = seen[static_cast<std::size_t>(old_index)];
        if (slot != 0U) {
            throw std::runtime_error("cell permutation is not bijective");
        }
        slot = 1U;
    }
}

void apply_cell_permutation(
    MeshTopology& topology,
    std::vector<Real>& volumes,
    const CellPermutation& permutation) {
    const Index n = static_cast<Index>(volumes.size());
    if (permutation.size() != n) {
        throw std::runtime_error("cell permutation size does not match volume count");
    }
    if (topology.coord_x.size() != static_cast<std::size_t>(n)
        || topology.coord_y.size() != static_cast<std::size_t>(n)
        || topology.is_boundary_cell.size() != static_cast<std::size_t>(n)
        || topology.neighbor_start.size() != static_cast<std::size_t>(n + 1)) {
        throw std::runtime_error("cell permutation cannot be applied to inconsistent topology");
    }

    MeshTopology out;
    out.coord_x.resize(static_cast<std::size_t>(n));
    out.coord_y.resize(static_cast<std::size_t>(n));
    out.is_boundary_cell.resize(static_cast<std::size_t>(n));
    std::vector<Real> reordered_volumes(static_cast<std::size_t>(n));

    for (Index new_index = 0; new_index < n; ++new_index) {
        const Index old_index = permutation.new_to_old()[static_cast<std::size_t>(new_index)];
        const auto ns = static_cast<std::size_t>(new_index);
        const auto os = static_cast<std::size_t>(old_index);
        out.coord_x[ns] = topology.coord_x[os];
        out.coord_y[ns] = topology.coord_y[os];
        out.is_boundary_cell[ns] = topology.is_boundary_cell[os];
        reordered_volumes[ns] = volumes[os];
    }

    out.neighbor_start.reserve(static_cast<std::size_t>(n + 1));
    out.neighbor_start.push_back(0);
    for (Index new_index = 0; new_index < n; ++new_index) {
        const Index old_index = permutation.new_to_old()[static_cast<std::size_t>(new_index)];
        const Index begin = topology.neighbor_start[static_cast<std::size_t>(old_index)];
        const Index end = topology.neighbor_start[static_cast<std::size_t>(old_index + 1)];
        for (Index edge = begin; edge < end; ++edge) {
            const auto es = static_cast<std::size_t>(edge);
            const Index old_neighbour = topology.neighbor_list[es];
            if (old_neighbour < 0 || old_neighbour >= n) {
                throw std::runtime_error("cell topology contains an out-of-range neighbour");
            }
            out.neighbor_list.push_back(permutation.old_to_new()[static_cast<std::size_t>(old_neighbour)]);
            out.interface_area.push_back(topology.interface_area[es]);
            out.distance.push_back(topology.distance[es]);
        }
        out.neighbor_start.push_back(static_cast<Index>(out.neighbor_list.size()));
    }

    topology = std::move(out);
    volumes = std::move(reordered_volumes);
}

} // namespace gss
