#pragma once

#include <span>
#include <vector>

#include "gss/core/Types.hpp"

namespace gss {

/// @brief Explicit cell permutation represented as new index to old index.
class CellPermutation {
public:
    CellPermutation() = default;

    /// @brief Build the identity permutation of size n.
    [[nodiscard]] static CellPermutation identity(Index n);

    /// @brief Build and validate a permutation stored as new index to old index.
    [[nodiscard]] static CellPermutation from_new_to_old(std::vector<Index> new_to_old);

    /// @brief Number of permuted cells.
    [[nodiscard]] Index size() const noexcept { return static_cast<Index>(new_to_old_.size()); }

    /// @brief Mapping from new index to old index.
    [[nodiscard]] std::span<const Index> new_to_old() const noexcept { return new_to_old_; }

    /// @brief Mapping from old index to new index.
    [[nodiscard]] std::span<const Index> old_to_new() const noexcept { return old_to_new_; }

    /// @brief True when the permutation leaves all indices unchanged.
    [[nodiscard]] bool is_identity() const noexcept;

    /// @brief Inverse permutation.
    [[nodiscard]] CellPermutation inverse() const;

    /// @brief Compose this permutation followed by a permutation of the new layout.
    [[nodiscard]] CellPermutation then(const CellPermutation& next) const;

private:
    explicit CellPermutation(std::vector<Index> new_to_old, std::vector<Index> old_to_new);

    std::vector<Index> new_to_old_{};
    std::vector<Index> old_to_new_{};
};

/// @brief Throw if values are not a bijection of [0,n).
void validate_cell_permutation(std::span<const Index> new_to_old);

/// @brief Apply a cell permutation to topology and cell volumes.
void apply_cell_permutation(
    MeshTopology& topology,
    std::vector<Real>& volumes,
    const CellPermutation& permutation);

} // namespace gss
