#pragma once

#include <span>

#include "gsmaker/core/Types.hpp"

namespace gss {

/**
 * @brief Finite-volume error norms for cell-centred scalar fields.
 *
 * The L1 and L2 norms are weighted by the cell volumes:
 *
 *   L1   = sum_i V_i |e_i|
 *   L2   = sqrt(sum_i V_i e_i^2)
 *   Linf = max_i |e_i|
 *
 * The volume-weighted RMS is
 *
 *   RMS = sqrt(sum_i V_i e_i^2 / sum_i V_i).
 *
 * No unweighted L2 overload is provided deliberately.  GSMakerLib supports
 * non-uniform unstructured meshes, for which treating every cell as having
 * equal measure would give the wrong discrete approximation to the
 * continuous L2 norm.
 */
class ErrorNorms {
public:
    struct Values {
        Real l1{};
        Real l2{};
        Real linf{};
        Real rms{};
        Real total_volume{};
    };

    /**
     * @brief Compute all supported norms in one pass.
     *
     * @throws std::invalid_argument if the arrays are empty, have different
     * sizes, contain non-finite values, or contain a non-positive cell volume.
     */
    [[nodiscard]] static Values compute(
        std::span<const Real> error,
        std::span<const Real> cell_volume);

    /**
     * @brief Volume-weighted discrete L1 norm.
     */
    [[nodiscard]] static Real l1(
        std::span<const Real> error,
        std::span<const Real> cell_volume);

    /**
     * @brief Volume-weighted discrete L2 norm.
     */
    [[nodiscard]] static Real l2(
        std::span<const Real> error,
        std::span<const Real> cell_volume);

    /**
     * @brief Maximum absolute error.
     *
     * Linf does not require volume weighting.
     *
     * @throws std::invalid_argument if the error array is empty or contains
     * non-finite values.
     */
    [[nodiscard]] static Real linf(
        std::span<const Real> error);

    /**
     * @brief Volume-weighted root-mean-square error.
     */
    [[nodiscard]] static Real rms(
        std::span<const Real> error,
        std::span<const Real> cell_volume);
};

} // namespace gss
