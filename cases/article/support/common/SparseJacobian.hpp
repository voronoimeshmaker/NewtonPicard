#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/core/Types.hpp"

namespace study {

using Real = gss::Real;
using Index = gss::Index;

/// @brief Sparse, graph-coloured assembly of B = d s_h / d Psi and J = A - B.
///
/// WHY THIS EXISTS.  The first implementation formed B by perturbing one cell
/// at a time, which costs 2n full source evaluations, and stored the result
/// densely, which costs n^2 entries and an O(n^3) factorisation.  At n = 768
/// that is 3072 whole-mesh gradient reconstructions per Jacobian and 589824
/// stored entries, and it is the reason every nonlinear campaign ran on 46 to
/// 855 cells while the frozen-Poisson study ran on 21760 to 100765.
///
/// The measured sparsity makes the cost unnecessary.  The stencil-radius study
/// reports dmax(B) = 1 and no entry beyond graph distance 2, on three mesh
/// families, two gradient methods and three resolutions: B has the pattern of
/// the face-neighbour graph.  Columns whose stencils do not overlap can then be
/// perturbed together, so the assembly costs 2 * (number of colours) source
/// evaluations rather than 2n, and the matrix holds a few entries per row.
///
/// THE COLOURING CONDITION.  Row i of B has entries at the closed r-ring of i.
/// Columns j and k may share a colour only if no row contains both, that is
/// only if their r-rings are disjoint, that is only if dist(j,k) > 2r.  The
/// conflict graph is therefore the mesh graph at distance 2r, and r is left as
/// a parameter so that the distance-1 measurement can be relied on by default
/// and widened to 2 whenever a new reconstruction operator is introduced.
///
/// VERIFY BEFORE TRUSTING.  assemble_dense_reference() keeps the original
/// one-column-at-a-time construction, and compare_to_dense() measures the
/// difference.  On any mesh small enough for both, they must agree to rounding;
/// a campaign that changes the gradient operator should run that comparison
/// once before relying on the sparse path.
namespace sparse {

/// Central-difference step for one component.  This mirrors
/// study::component_central_step in FadeevNonlinearStudyCommon.hpp and must
/// stay identical to it, so that the sparse and dense assemblies differ only
/// in how the perturbations are grouped.
[[nodiscard]] inline Real central_step(
    Real psi_value,
    Real scale = Real{1}) {

    const Real epsilon =
        static_cast<Real>(std::numeric_limits<double>::epsilon());

    const Real characteristic =
        static_cast<Real>(std::cbrt(static_cast<double>(epsilon)));

    const Real magnitude =
        static_cast<Real>(std::abs(static_cast<double>(psi_value)));

    return scale * characteristic * std::max(Real{1}, magnitude);
}

/// Closed neighbourhood of every cell out to graph distance `radius`,
/// sorted ascending, diagonal included.
[[nodiscard]] inline std::vector<std::vector<Index>> closed_rings(
    const gss::MeshTopology& topology,
    std::size_t radius) {

    const std::size_t n = topology.coord_x.size();

    if (n == 0U) {
        throw std::invalid_argument("closed_rings requires a non-empty mesh");
    }
    if (radius == 0U) {
        throw std::invalid_argument("closed_rings requires a positive radius");
    }

    std::vector<std::vector<Index>> ring(n);
    std::vector<std::size_t> visited(n, n + 1U);
    std::vector<Index> frontier;
    std::vector<Index> next;

    for (std::size_t seed = 0; seed < n; ++seed) {

        frontier.clear();
        frontier.push_back(static_cast<Index>(seed));
        visited[seed] = seed;

        std::vector<Index>& target = ring[seed];
        target.clear();
        target.push_back(static_cast<Index>(seed));

        for (std::size_t level = 0; level < radius; ++level) {

            next.clear();

            for (const Index current : frontier) {

                const auto c = static_cast<std::size_t>(current);
                const auto begin = static_cast<std::size_t>(topology.neighbor_start[c]);
                const auto end = static_cast<std::size_t>(topology.neighbor_start[c + 1U]);

                for (std::size_t k = begin; k < end; ++k) {

                    const auto neighbour =
                        static_cast<std::size_t>(topology.neighbor_list[k]);

                    if (neighbour >= n) {
                        continue;
                    }
                    if (visited[neighbour] == seed) {
                        continue;
                    }

                    visited[neighbour] = seed;
                    target.push_back(static_cast<Index>(neighbour));
                    next.push_back(static_cast<Index>(neighbour));
                }
            }

            frontier.swap(next);
        }

        std::sort(target.begin(), target.end());
    }

    return ring;
}

/// Sparsity pattern of B together with a colouring valid for
/// central-difference assembly.
struct Pattern {
    std::size_t cells{};
    std::size_t radius{};

    /// CSR structure of B: row_start has cells + 1 entries.
    std::vector<Index> row_start;
    std::vector<Index> column;

    /// Colour of every cell, in [0, colours).
    std::vector<Index> colour;
    std::size_t colours{};

    /// Cells belonging to each colour.
    std::vector<std::vector<Index>> colour_members;

    [[nodiscard]] std::size_t nnz() const noexcept { return column.size(); }
};

/// Build the pattern and colour it.
///
/// The colouring is greedy over cells ordered by decreasing conflict degree,
/// which is the standard heuristic and gives a colour count close to the
/// maximum degree plus one.  Any valid colouring gives a correct Jacobian; the
/// ordering only affects how many source evaluations the assembly costs.
[[nodiscard]] inline Pattern build_pattern(
    const gss::MeshTopology& topology,
    std::size_t radius = 1U) {

    const std::size_t n = topology.coord_x.size();

    Pattern pattern;
    pattern.cells = n;
    pattern.radius = radius;

    const auto ring = closed_rings(topology, radius);

    pattern.row_start.assign(n + 1U, Index{0});
    pattern.column.clear();

    for (std::size_t i = 0; i < n; ++i) {
        for (const Index j : ring[i]) {
            pattern.column.push_back(j);
        }
        pattern.row_start[i + 1U] = static_cast<Index>(pattern.column.size());
    }

    // Columns conflict when their rings overlap, that is at mesh distance
    // 2 * radius or less.
    const auto conflict = closed_rings(topology, 2U * radius);

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});

    std::sort(
        order.begin(),
        order.end(),
        [&conflict](std::size_t a, std::size_t b) {
            if (conflict[a].size() != conflict[b].size()) {
                return conflict[a].size() > conflict[b].size();
            }
            return a < b;
        });

    constexpr Index uncoloured = -1;

    pattern.colour.assign(n, uncoloured);

    std::vector<unsigned char> taken;

    for (const std::size_t cell : order) {

        taken.assign(conflict[cell].size() + 1U, 0U);

        for (const Index other : conflict[cell]) {

            const auto o = static_cast<std::size_t>(other);

            if (o == cell) {
                continue;
            }
            if (pattern.colour[o] == uncoloured) {
                continue;
            }

            const auto used = static_cast<std::size_t>(pattern.colour[o]);

            if (used < taken.size()) {
                taken[used] = 1U;
            }
        }

        Index chosen = 0;

        while (static_cast<std::size_t>(chosen) < taken.size()) {
            if (taken[static_cast<std::size_t>(chosen)] == 0U) {
                break;
            }
            ++chosen;
        }

        pattern.colour[cell] = chosen;
        pattern.colours = std::max(
            pattern.colours,
            static_cast<std::size_t>(chosen) + 1U);
    }

    pattern.colour_members.assign(pattern.colours, {});

    for (std::size_t i = 0; i < n; ++i) {
        const auto c = static_cast<std::size_t>(pattern.colour[i]);
        pattern.colour_members[c].push_back(static_cast<Index>(i));
    }

    return pattern;
}

/// Assemble B in the given pattern by coloured central differences.
///
/// `source` maps a state to the integrated nonlinear source, exactly as
/// study::evaluate_integrated_source does.  It is called 2 * colours times.
template<class SourceEvaluator>
[[nodiscard]] std::vector<Real> assemble_values(
    const Pattern& pattern,
    std::span<const Real> state,
    SourceEvaluator&& source,
    Real finite_difference_scale = Real{1}) {

    const std::size_t n = pattern.cells;

    if (state.size() != n) {
        throw std::invalid_argument("assemble_values state size does not match the pattern");
    }

    std::vector<Real> values(pattern.nnz(), Real{0});

    std::vector<Real> step(n, Real{0});

    for (std::size_t i = 0; i < n; ++i) {
        step[i] = central_step(state[i], finite_difference_scale);
    }

    std::vector<Real> plus(state.begin(), state.end());
    std::vector<Real> minus(state.begin(), state.end());

    for (std::size_t c = 0; c < pattern.colours; ++c) {

        const std::vector<Index>& members = pattern.colour_members[c];

        for (const Index member : members) {
            const auto j = static_cast<std::size_t>(member);
            plus[j] = state[j] + step[j];
            minus[j] = state[j] - step[j];
        }

        const auto source_plus = source(std::span<const Real>(plus));
        const auto source_minus = source(std::span<const Real>(minus));

        if (source_plus.size() != n || source_minus.size() != n) {
            throw std::runtime_error("assemble_values source evaluator returned the wrong size");
        }

        // Exactly one member of this colour can influence any given row, by
        // construction of the colouring, so the difference is attributed
        // unambiguously.
        for (std::size_t i = 0; i < n; ++i) {

            const auto begin = static_cast<std::size_t>(pattern.row_start[i]);
            const auto end = static_cast<std::size_t>(pattern.row_start[i + 1U]);

            for (std::size_t k = begin; k < end; ++k) {

                const auto j = static_cast<std::size_t>(pattern.column[k]);

                if (static_cast<std::size_t>(pattern.colour[j]) != c) {
                    continue;
                }

                const Real denominator = Real{2} * step[j];

                values[k] = (source_plus[i] - source_minus[i]) / denominator;
            }
        }

        for (const Index member : members) {
            const auto j = static_cast<std::size_t>(member);
            plus[j] = state[j];
            minus[j] = state[j];
        }
    }

    return values;
}

/// B as a CsrMatrix in the pattern.
[[nodiscard]] inline gss::CsrMatrix to_csr(
    const Pattern& pattern,
    std::vector<Real> values) {

    if (values.size() != pattern.nnz()) {
        throw std::invalid_argument("to_csr value count does not match the pattern");
    }

    return gss::CsrMatrix{
        static_cast<Index>(pattern.cells),
        pattern.row_start,
        pattern.column,
        std::move(values)};
}

/// J = A - B, assembled in the union of the two patterns.
///
/// A is the finite-volume Poisson operator, whose pattern is the face-neighbour
/// graph; for radius 1 or greater the pattern of B already contains it, but the
/// union is formed explicitly so that the routine stays correct if either
/// operator changes.
[[nodiscard]] inline gss::CsrMatrix subtract_from(
    const gss::CsrMatrix& a,
    const Pattern& pattern,
    std::span<const Real> b_values) {

    const auto n = static_cast<std::size_t>(a.rows());

    if (n != pattern.cells) {
        throw std::invalid_argument("subtract_from matrix and pattern sizes differ");
    }
    if (b_values.size() != pattern.nnz()) {
        throw std::invalid_argument("subtract_from value count does not match the pattern");
    }

    const auto a_row = a.row_start();
    const auto a_col = a.column();
    const auto a_val = a.value();

    std::vector<Index> row_start(n + 1U, Index{0});
    std::vector<Index> column;
    std::vector<Real> value;

    column.reserve(pattern.nnz() + a.nnz());
    value.reserve(pattern.nnz() + a.nnz());

    std::vector<Index> merged;

    for (std::size_t i = 0; i < n; ++i) {

        merged.clear();

        const auto ab = static_cast<std::size_t>(a_row[i]);
        const auto ae = static_cast<std::size_t>(a_row[i + 1U]);
        const auto bb = static_cast<std::size_t>(pattern.row_start[i]);
        const auto be = static_cast<std::size_t>(pattern.row_start[i + 1U]);

        for (std::size_t k = ab; k < ae; ++k) {
            merged.push_back(a_col[k]);
        }
        for (std::size_t k = bb; k < be; ++k) {
            merged.push_back(pattern.column[k]);
        }

        std::sort(merged.begin(), merged.end());
        merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

        for (const Index j : merged) {

            Real entry = Real{0};

            for (std::size_t k = ab; k < ae; ++k) {
                if (a_col[k] == j) {
                    entry += a_val[k];
                    break;
                }
            }
            for (std::size_t k = bb; k < be; ++k) {
                if (pattern.column[k] == j) {
                    entry -= b_values[k];
                    break;
                }
            }

            column.push_back(j);
            value.push_back(entry);
        }

        row_start[i + 1U] = static_cast<Index>(value.size());
    }

    return gss::CsrMatrix{
        static_cast<Index>(n),
        std::move(row_start),
        std::move(column),
        std::move(value)};
}

/// The original one-column-at-a-time construction, kept as the reference the
/// sparse path is verified against.  Cost is 2n source evaluations and n^2
/// storage, so this is for small meshes only.
template<class SourceEvaluator>
[[nodiscard]] std::vector<Real> assemble_dense_reference(
    std::span<const Real> state,
    SourceEvaluator&& source,
    Real finite_difference_scale = Real{1}) {

    const std::size_t n = state.size();

    std::vector<Real> dense(n * n, Real{0});
    std::vector<Real> plus(state.begin(), state.end());
    std::vector<Real> minus(state.begin(), state.end());

    for (std::size_t j = 0; j < n; ++j) {

        const Real h = central_step(state[j], finite_difference_scale);

        plus[j] = state[j] + h;
        minus[j] = state[j] - h;

        const auto source_plus = source(std::span<const Real>(plus));
        const auto source_minus = source(std::span<const Real>(minus));

        const Real denominator = Real{2} * h;

        for (std::size_t i = 0; i < n; ++i) {
            dense[i * n + j] = (source_plus[i] - source_minus[i]) / denominator;
        }

        plus[j] = state[j];
        minus[j] = state[j];
    }

    return dense;
}

struct PatternComparison {
    /// Largest difference on entries the sparse pattern keeps.
    double max_difference_inside{0.0};

    /// Largest magnitude the sparse pattern discards.  A pattern that is too
    /// narrow shows up here and nowhere else.
    double max_magnitude_outside{0.0};

    /// Largest magnitude anywhere in the dense reference.
    double max_magnitude{0.0};

    [[nodiscard]] double relative_inside() const {
        return max_magnitude > 0.0 ? max_difference_inside / max_magnitude : 0.0;
    }

    [[nodiscard]] double relative_outside() const {
        return max_magnitude > 0.0 ? max_magnitude_outside / max_magnitude : 0.0;
    }
};

/// Compare a coloured assembly against the dense reference.
///
/// relative_inside measures whether the colouring attributed the differences
/// correctly; relative_outside measures whether the pattern is wide enough.
/// The second is the one that matters when a new gradient operator is
/// introduced, because a narrow pattern is silently wrong rather than noisy.
[[nodiscard]] inline PatternComparison compare_to_dense(
    const Pattern& pattern,
    std::span<const Real> sparse_values,
    const std::vector<Real>& dense) {

    const std::size_t n = pattern.cells;

    if (dense.size() != n * n) {
        throw std::invalid_argument("compare_to_dense reference has the wrong size");
    }
    if (sparse_values.size() != pattern.nnz()) {
        throw std::invalid_argument("compare_to_dense value count does not match the pattern");
    }

    PatternComparison comparison;

    std::vector<unsigned char> inside(n, 0U);

    for (std::size_t i = 0; i < n; ++i) {

        std::fill(inside.begin(), inside.end(), 0U);

        const auto begin = static_cast<std::size_t>(pattern.row_start[i]);
        const auto end = static_cast<std::size_t>(pattern.row_start[i + 1U]);

        for (std::size_t k = begin; k < end; ++k) {

            const auto j = static_cast<std::size_t>(pattern.column[k]);

            inside[j] = 1U;

            const double reference = static_cast<double>(dense[i * n + j]);
            const double got = static_cast<double>(sparse_values[k]);

            comparison.max_difference_inside = std::max(
                comparison.max_difference_inside,
                std::abs(got - reference));
        }

        for (std::size_t j = 0; j < n; ++j) {

            const double reference = std::abs(static_cast<double>(dense[i * n + j]));

            comparison.max_magnitude = std::max(comparison.max_magnitude, reference);

            if (inside[j] == 0U) {
                comparison.max_magnitude_outside = std::max(
                    comparison.max_magnitude_outside,
                    reference);
            }
        }
    }

    return comparison;
}

} // namespace sparse
} // namespace study