#pragma once

#include <cmath>
#include <stdexcept>

#include "gsmaker/core/Config.hpp"
#include "gsmaker/core/Types.hpp"

namespace gss {

/// @brief Configuration for relaxed Picard fixed-point iteration.
///
/// This object contains nonlinear iteration parameters only. Linear KSP/PC
/// configuration remains owned by PetscLinearConfig.
struct PicardConfig {
    Index max_iterations{100};

    Real relative_residual_tolerance{Real{1.0e-8}};
    Real absolute_residual_tolerance{Real{1.0e-12}};

    /// Relaxation parameter in
    /// psi_{k+1} = (1-omega) psi_k + omega psi_hat_{k+1}.
    /// Picard itself imposes only omega > 0; individual studies may restrict
    /// the admissible range further.
    Real omega{Real{1}};

    [[nodiscard]]
    static PicardConfig from_config(const Config& config) {
        PicardConfig result;
        result.max_iterations = config.solver.max_iter;
        result.relative_residual_tolerance = config.solver.tolerance_relative;
        result.absolute_residual_tolerance = config.solver.tolerance_absolute;
        result.omega = config.solver.relaxation_omega;
        result.validate();
        return result;
    }

    void validate() const {
        if (max_iterations <= 0) {
            throw std::invalid_argument("PicardConfig max_iterations must be positive");
        }
        if (!std::isfinite(static_cast<double>(relative_residual_tolerance)) ||
            relative_residual_tolerance < Real{0}) {
            throw std::invalid_argument(
                "PicardConfig relative residual tolerance must be finite and non-negative");
        }
        if (!std::isfinite(static_cast<double>(absolute_residual_tolerance)) ||
            absolute_residual_tolerance < Real{0}) {
            throw std::invalid_argument(
                "PicardConfig absolute residual tolerance must be finite and non-negative");
        }
        if (!std::isfinite(static_cast<double>(omega)) || !(omega > Real{0})) {
            throw std::invalid_argument("PicardConfig omega must be finite and positive");
        }
    }
};

} // namespace gss
