#pragma once

#include <cstddef>
#include <vector>

#include "gsmaker/core/Types.hpp"


namespace gss {

/// @brief Sparsity pattern of a stencil operator together with a colouring
/// valid for finite-difference assembly.
///
/// The pattern is stored in CSR structure so that it can be handed directly to
/// an assembler or to CsrMatrix.  Columns carrying the same colour may be
/// perturbed simultaneously without any row receiving two contributions.
struct SparsityPattern {
    Index rows{0};

    /// Graph distance used to build the stencil of each row.
    Index stencil_radius{0};

    std::vector<Index> row_start{};
    std::vector<Index> column{};

    /// Colour of every row index, in [0, colours).
    std::vector<Index> colour{};
    Index colours{0};

    /// Row indices belonging to each colour.
    std::vector<std::vector<Index>> colour_members{};

    [[nodiscard]]
    std::size_t nnz() const noexcept {
        return column.size();
    }
};


/// @brief Build and colour the sparsity pattern of a mesh stencil operator.
///
/// An operator whose row i depends only on the cells within graph distance r of
/// i can be assembled by finite differences in far fewer than n perturbations.
/// Columns j and k may be perturbed together only if no row contains both, that
/// is only if their r-rings are disjoint, that is only if dist(j,k) exceeds 2r.
/// The conflict graph is therefore the mesh graph at distance 2r, and the
/// number of colours is bounded by its maximum degree rather than by the cell
/// count.
///
/// The pattern depends only on the mesh, so it is preparation that belongs
/// before any nonlinear iteration and is reused unchanged by every assembly.
class SparsityAnalysis final {
public:
    /// Closed neighbourhood of every cell out to the given graph distance,
    /// sorted ascending, with the cell itself included.
    [[nodiscard]]
    static std::vector<std::vector<Index>> closed_rings(
        const MeshTopology& topology,
        Index radius);

    /// Pattern of an operator with the given stencil radius, together with a
    /// colouring valid for central-difference assembly.
    ///
    /// The colouring is greedy over rows ordered by decreasing conflict degree.
    /// Any valid colouring yields a correct assembly; the ordering affects only
    /// how many colours, and therefore how many evaluations, are required.
    [[nodiscard]]
    static SparsityPattern build(
        const MeshTopology& topology,
        Index radius = Index{1});

    /// Verify that no two rows sharing a colour lie within distance 2r of one
    /// another.  A colouring that fails this is silently wrong rather than
    /// noisy, so it is checked rather than assumed.
    [[nodiscard]]
    static bool colouring_is_valid(
        const MeshTopology& topology,
        const SparsityPattern& pattern);
};

} // namespace gss