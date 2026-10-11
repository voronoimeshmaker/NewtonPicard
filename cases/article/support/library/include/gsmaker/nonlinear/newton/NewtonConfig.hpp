#pragma once

#include <cmath>
#include <functional>
#include <span>
#include <stdexcept>

#include "gsmaker/core/Types.hpp"

namespace gss {

/// @brief Numerical controls for the generic damped Newton algorithm.
///
/// This type is intentionally distinct from gss::NewtonConfig in
/// gsmaker/core/Config.hpp. The latter belongs to the global/YAML problem
/// configuration and currently stores problem-level Newton options such as
/// finite-difference Jacobian and line-search settings. NewtonSolverConfig
/// contains only the controls consumed directly by gss::Newton::solve().
struct NewtonSolverConfig {
    Index max_iterations{50};

    Real relative_residual_tolerance{Real{1.0e-8}};
    Real absolute_residual_tolerance{Real{1.0e-12}};

    /// Constant damping applied after solving J delta = -R:
    /// psi_{k+1} = psi_k + damping * delta.
    Real damping{Real{1}};

    /// Safeguard, NOT a globalisation.
    ///
    /// A step is shortened only when the resulting state cannot be used at
    /// all: a non-finite component, a state the admissibility predicate
    /// rejects, or a residual evaluation that fails or returns a non-finite
    /// value. A step that merely increases the residual is accepted, because
    /// requiring decrease would make this a line search and change the method
    /// under study.
    Index max_step_reductions{10};

    /// Factor applied to the step length at each reduction.
    Real step_reduction_factor{Real{0.5}};

    /// Optional domain test, evaluated on each trial state. Empty means every
    /// finite state is admissible.
    ///
    /// For the Grad-Shafranov source this is where the analytical bound
    /// belongs: Psi is bounded below by Fadeev::psi_lower_bound(), and an
    /// iterate far beneath it has left the set on which the problem is posed
    /// several hundred iterations before the arithmetic degrades.
    std::function<bool(std::span<const Real>)> admissible{};

    void validate() const {
        if (max_iterations <= 0) {
            throw std::invalid_argument(
                "NewtonSolverConfig max_iterations must be positive");
        }
        if (!std::isfinite(static_cast<double>(relative_residual_tolerance)) ||
            relative_residual_tolerance < Real{0}) {
            throw std::invalid_argument(
                "NewtonSolverConfig relative residual tolerance must be finite and non-negative");
        }
        if (!std::isfinite(static_cast<double>(absolute_residual_tolerance)) ||
            absolute_residual_tolerance < Real{0}) {
            throw std::invalid_argument(
                "NewtonSolverConfig absolute residual tolerance must be finite and non-negative");
        }
        if (!std::isfinite(static_cast<double>(damping)) ||
            !(damping > Real{0}) || damping > Real{1}) {
            throw std::invalid_argument(
                "NewtonSolverConfig damping must be finite and satisfy 0 < damping <= 1");
        }
        if (max_step_reductions < 0) {
            throw std::invalid_argument(
                "NewtonSolverConfig max_step_reductions must be non-negative");
        }
        if (!std::isfinite(static_cast<double>(step_reduction_factor)) ||
            !(step_reduction_factor > Real{0}) ||
            step_reduction_factor >= Real{1}) {
            throw std::invalid_argument(
                "NewtonSolverConfig step reduction factor must satisfy 0 < factor < 1");
        }
    }
};

} // namespace gss