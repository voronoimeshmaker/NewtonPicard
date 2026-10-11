#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "gss/core/Types.hpp"

namespace gss {

/// @brief Bundle of area-weighted L1/L2 and unweighted Linf norms.
struct NormBundle {
    Real L1{0};
    Real L2{0};
    Real Linf{0};
    Real V_total{0};
};

/// @brief Compute area-weighted L1/L2 and unweighted Linf norms.
[[nodiscard]] inline NormBundle compute_norms(
    std::span<const Real> error_field,
    std::span<const Real> cell_volumes) noexcept {
    NormBundle out{};
    const auto n = std::min(error_field.size(), cell_volumes.size());
    for (std::size_t i = 0; i < n; ++i) {
        const Real e = std::abs(error_field[i]);
        const Real v = cell_volumes[i];
        out.L1 += e * v;
        out.L2 += e * e * v;
        out.Linf = std::max(out.Linf, e);
        out.V_total += v;
    }
    if (out.V_total > 0) {
        out.L1 /= out.V_total;
        out.L2 = std::sqrt(out.L2 / out.V_total);
    }
    return out;
}

/// @brief Norms split into bulk and island cells.
struct SeparatedNorms {
    NormBundle bulk;
    NormBundle island;
    Index n_bulk_cells{0};
    Index n_island_cells{0};
};

/// @brief Compute norms separately for cells below and above a kappa threshold.
[[nodiscard]] inline SeparatedNorms compute_separated_norms(
    std::span<const Real> error_field,
    std::span<const Real> cell_volumes,
    std::span<const Real> kappa,
    Real kappa_island_threshold) noexcept {
    SeparatedNorms out{};
    std::vector<Real> bulk_e;
    std::vector<Real> bulk_v;
    std::vector<Real> island_e;
    std::vector<Real> island_v;
    const auto n = std::min({error_field.size(), cell_volumes.size(), kappa.size()});
    bulk_e.reserve(n);
    bulk_v.reserve(n);
    island_e.reserve(n);
    island_v.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (kappa[i] > kappa_island_threshold) {
            island_e.push_back(error_field[i]);
            island_v.push_back(cell_volumes[i]);
            ++out.n_island_cells;
        } else {
            bulk_e.push_back(error_field[i]);
            bulk_v.push_back(cell_volumes[i]);
            ++out.n_bulk_cells;
        }
    }
    out.bulk = compute_norms(bulk_e, bulk_v);
    out.island = compute_norms(island_e, island_v);
    return out;
}

} // namespace gss
