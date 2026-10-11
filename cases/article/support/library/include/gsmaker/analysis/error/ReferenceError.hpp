#pragma once

#include <span>
#include <vector>

#include "gsmaker/core/Types.hpp"

namespace gss {

/**
 * @brief Construct the pointwise error between a numerical field and a
 * reference field.
 *
 * The error convention is
 *
 *   e_i = numerical_i - reference_i.
 *
 * This class deliberately does not compute norms.  Norm evaluation belongs to
 * ErrorNorms, which applies the appropriate finite-volume weights.
 */
class ReferenceError {
public:
    /**
     * @brief Compute the pointwise reference error.
     *
     * @param numerical Cell-centred numerical values.
     * @param reference Cell-centred reference values at the same locations.
     *
     * @return Vector containing numerical - reference for every entry.
     *
     * @throws std::invalid_argument if the arrays are empty, have different
     * sizes, or contain non-finite values.
     */
    [[nodiscard]] static std::vector<Real> compute(
        std::span<const Real> numerical,
        std::span<const Real> reference);
};

} // namespace gss
