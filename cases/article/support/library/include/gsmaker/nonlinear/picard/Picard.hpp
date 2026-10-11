#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <limits>
#include <functional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/core/Types.hpp"
#include "gsmaker/linear_system/PetscLinearSystemController.hpp"
#include "gsmaker/nonlinear/picard/PicardConfig.hpp"
#include "gsmaker/nonlinear/picard/PicardResult.hpp"


namespace gss {

/// @brief Per-iterate data passed synchronously to an optional observer.
///
/// state is valid only for the duration of the observer call. Iteration zero
/// represents the supplied initial state and has zero step/linear statistics.
struct PicardIterationInfo {
    Index iteration{0};
    std::span<const Real> state{};

    Real residual_norm{0};
    Real relative_residual{0};
    Real step_norm{0};
    Real relative_step{0};

    Index linear_iterations{0};
    Real linear_true_relative_residual{0};
};


struct NoPicardObserver {
    void operator()(const PicardIterationInfo&) const noexcept {}
};

namespace detail {

[[nodiscard]]
inline Real picard_norm2(std::span<const Real> values) {
    long double sum = 0.0L;
    for (const Real value : values) {
        const long double v = static_cast<long double>(value);
        sum += v * v;
    }
    return static_cast<Real>(std::sqrt(sum));
}

[[nodiscard]]
inline Real picard_relative_denominator(Real value) {
    const Real floor = static_cast<Real>(
        100.0 * std::numeric_limits<double>::epsilon());
    return std::max(value, floor);
}

[[nodiscard]]
inline std::vector<Real> picard_residual(
    const CsrMatrix& matrix,
    std::span<const Real> state,
    std::span<const Real> nonlinear_rhs,
    std::span<const Real> boundary_rhs) {

    if (nonlinear_rhs.size() != state.size() ||
        boundary_rhs.size() != state.size()) {
        throw std::invalid_argument("Picard residual vector size mismatch");
    }

    auto residual = matrix.multiply(state);
    for (std::size_t i = 0; i < residual.size(); ++i) {
        residual[i] -= nonlinear_rhs[i] + boundary_rhs[i];
    }
    return residual;
}

template<class SourceEvaluator>
[[nodiscard]]
std::vector<Real> evaluate_picard_source(
    SourceEvaluator& source,
    std::span<const Real> state) {

    auto values = source(state);
    static_assert(
        std::same_as<std::remove_cvref_t<decltype(values)>, std::vector<Real>>,
        "Picard source evaluator must return std::vector<gss::Real>");

    if (values.size() != state.size()) {
        throw std::invalid_argument("Picard source evaluator returned the wrong vector size");
    }
    for (const Real value : values) {
        if (!std::isfinite(static_cast<double>(value))) {
            throw std::runtime_error("Picard source evaluator returned a non-finite value");
        }
    }
    return values;
}

} // namespace detail

/// @brief Generic relaxed Picard iteration for A psi = s(psi) + b.
///
/// Picard does not know how the nonlinear source is evaluated. The caller
/// supplies a callable state -> integrated nonlinear RHS. The matrix and the
/// PETSc controller are reused unchanged throughout the iteration.
class Picard final {
public:
    template<class SourceEvaluator, class Observer = NoPicardObserver>
    [[nodiscard]]
    static PicardResult solve(
        const CsrMatrix& matrix,
        std::span<const Real> boundary_rhs,
        std::span<const Real> initial_state,
        SourceEvaluator&& source_evaluator,
        PetscLinearSystemController& linear_solver,
        PicardConfig config,
        Observer&& observer = {}) {

        using clock = std::chrono::steady_clock;
        const auto start = clock::now();

        config.validate();
        matrix.validate();

        const auto n = static_cast<std::size_t>(matrix.rows());
        if (n == 0U || initial_state.size() != n || boundary_rhs.size() != n) {
            throw std::invalid_argument("Picard matrix and vector dimensions must agree and be non-zero");
        }
        if (linear_solver.rows() != matrix.rows()) {
            throw std::invalid_argument("Picard linear solver dimension does not match matrix");
        }

        PicardResult result;
        result.solution.assign(initial_state.begin(), initial_state.end());

        auto source = detail::evaluate_picard_source(
            source_evaluator,
            std::span<const Real>(result.solution));
        auto residual = detail::picard_residual(
            matrix,
            result.solution,
            source,
            boundary_rhs);

        const Real initial_residual = detail::picard_norm2(residual);
        const Real residual_scale = detail::picard_relative_denominator(initial_residual);

        result.initial_residual_norm = initial_residual;
        result.final_residual_norm = initial_residual;
        result.final_relative_residual = initial_residual / residual_scale;
        result.residual_history.push_back(initial_residual);
        result.relative_residual_history.push_back(result.final_relative_residual);

        std::invoke(
            observer,
            PicardIterationInfo{
                0,
                std::span<const Real>(result.solution),
                initial_residual,
                result.final_relative_residual,
                Real{0},
                Real{0},
                0,
                Real{0}});

        if (initial_residual <= config.absolute_residual_tolerance) {
            result.converged = true;
            result.termination = PicardTermination::ConvergedAbsoluteResidual;
            result.total_seconds = std::chrono::duration<double>(clock::now() - start).count();
            return result;
        }

        std::vector<Real> rhs(n, Real{0});
        std::vector<Real> previous = result.solution;

        for (Index iteration = 0; iteration < config.max_iterations; ++iteration) {
            for (std::size_t i = 0; i < n; ++i) {
                rhs[i] = source[i] + boundary_rhs[i];
            }

            // Solve the frozen linear problem from the controller's default
            // zero initial guess.  In particular, PETSc KSPPREONLY rejects a
            // non-zero initial guess because it applies the preconditioner
            // exactly once rather than performing an iterative correction.
            // Reusing A/KSP/PC still preserves the important Picard setup
            // reuse; only the solution vector is reinitialised for each RHS.
            const auto linear = linear_solver.solve(rhs);
            ++result.linear_solves;
            result.total_linear_iterations += linear.iterations;
            result.max_linear_true_relative_residual = std::max(
                result.max_linear_true_relative_residual,
                static_cast<Real>(linear.true_relative_residual));

            if (!linear.converged()) {
                result.termination = PicardTermination::LinearSolveFailure;
                result.linear_failure_reason = linear.reason_string;
                result.total_seconds = std::chrono::duration<double>(clock::now() - start).count();
                return result;
            }

            previous = result.solution;
            for (std::size_t i = 0; i < n; ++i) {
                result.solution[i] =
                    (Real{1} - config.omega) * previous[i] +
                    config.omega * linear.solution[i];
            }

            std::vector<Real> step(n, Real{0});
            for (std::size_t i = 0; i < n; ++i) {
                step[i] = result.solution[i] - previous[i];
            }

            const Real step_norm = detail::picard_norm2(step);
            const Real state_norm = detail::picard_norm2(result.solution);
            const Real relative_step =
                step_norm / std::max(Real{1}, detail::picard_relative_denominator(state_norm));

            source = detail::evaluate_picard_source(
                source_evaluator,
                std::span<const Real>(result.solution));
            residual = detail::picard_residual(
                matrix,
                result.solution,
                source,
                boundary_rhs);

            const Real residual_norm = detail::picard_norm2(residual);
            const Real relative_residual = residual_norm / residual_scale;

            result.iterations = iteration + 1;
            result.final_residual_norm = residual_norm;
            result.final_relative_residual = relative_residual;
            result.final_step_norm = step_norm;
            result.final_relative_step = relative_step;

            result.residual_history.push_back(residual_norm);
            result.relative_residual_history.push_back(relative_residual);
            result.step_norm_history.push_back(step_norm);
            result.relative_step_history.push_back(relative_step);
            result.linear_iterations_history.push_back(linear.iterations);
            result.linear_true_relative_residual_history.push_back(
                static_cast<Real>(linear.true_relative_residual));

            std::invoke(
                observer,
                PicardIterationInfo{
                    iteration + 1,
                    std::span<const Real>(result.solution),
                    residual_norm,
                    relative_residual,
                    step_norm,
                    relative_step,
                    linear.iterations,
                    static_cast<Real>(linear.true_relative_residual)});

            if (residual_norm <= config.absolute_residual_tolerance) {
                result.converged = true;
                result.termination = PicardTermination::ConvergedAbsoluteResidual;
                break;
            }
            if (relative_residual <= config.relative_residual_tolerance) {
                result.converged = true;
                result.termination = PicardTermination::ConvergedRelativeResidual;
                break;
            }
        }

        result.total_seconds = std::chrono::duration<double>(clock::now() - start).count();
        return result;
    }
};

} // namespace gss
