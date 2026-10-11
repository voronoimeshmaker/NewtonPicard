#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/core/Types.hpp"
#include "gsmaker/linear_system/PetscLinearConfig.hpp"
#include "gsmaker/linear_system/PetscLinearSystemController.hpp"
#include "gsmaker/nonlinear/newton/NewtonSolverConfig.hpp"
#include "gsmaker/nonlinear/newton/NewtonResult.hpp"

namespace gss {

/// @brief Data supplied synchronously to an optional Newton observer.
///
/// iteration zero is the supplied initial state. The state span is valid only
/// during the observer call.
struct NewtonIterationInfo {
    Index iteration{0};
    std::span<const Real> state{};

    Real residual_norm{0};
    Real relative_residual{0};
    Real step_norm{0};
    Real relative_step{0};
    Real damping{0};

    Index linear_iterations{0};
    Real linear_true_relative_residual{0};

    Index step_reductions{0};
};

struct NoNewtonObserver {
    void operator()(const NewtonIterationInfo&) const noexcept {}
};

namespace detail {

[[nodiscard]]
inline Real newton_norm2(std::span<const Real> values) {
    long double sum = 0.0L;
    for (const Real value : values) {
        const long double v = static_cast<long double>(value);
        sum += v * v;
    }
    return static_cast<Real>(std::sqrt(sum));
}

[[nodiscard]]
inline Real newton_relative_denominator(Real value) {
    const Real floor = static_cast<Real>(
        100.0 * std::numeric_limits<double>::epsilon());
    return std::max(value, floor);
}

[[nodiscard]]
inline bool newton_all_finite(std::span<const Real> values) noexcept {
    for (const Real value : values) {
        if (!std::isfinite(static_cast<double>(value))) {
            return false;
        }
    }
    return true;
}

/// Outcome of attempting a residual evaluation that is permitted to fail.
struct NewtonResidualAttempt {
    std::vector<Real> values;
    bool ok{false};
    std::string reason;
};

/// Evaluate the residual without letting a failure end the solve.
///
/// The evaluator may fail in two ways that a nonlinear solver must survive:
/// it may return non-finite values, and it may throw, which is what the
/// gradient-reconstruction source does when the state it is handed is not
/// usable. Both are reported, never propagated.
template<class ResidualEvaluator>
[[nodiscard]]
NewtonResidualAttempt try_evaluate_newton_residual(
    ResidualEvaluator& residual_evaluator,
    std::span<const Real> state) {

    NewtonResidualAttempt attempt;

    try {
        auto values = residual_evaluator(state);
        static_assert(
            std::same_as<std::remove_cvref_t<decltype(values)>, std::vector<Real>>,
            "Newton residual evaluator must return std::vector<gss::Real>");

        if (values.size() != state.size()) {
            attempt.reason =
                "residual evaluator returned the wrong vector size";
            return attempt;
        }

        if (!newton_all_finite(values)) {
            attempt.reason =
                "residual evaluator returned a non-finite value";
            return attempt;
        }

        attempt.values = std::move(values);
        attempt.ok = true;
        return attempt;

    } catch (const std::exception& error) {
        attempt.reason = error.what();
        return attempt;
    }
}

template<class JacobianEvaluator>
[[nodiscard]]
std::optional<CsrMatrix> try_evaluate_newton_jacobian(
    JacobianEvaluator& jacobian_evaluator,
    std::span<const Real> state,
    std::string& reason) {

    try {
        auto jacobian = jacobian_evaluator(state);
        static_assert(
            std::same_as<std::remove_cvref_t<decltype(jacobian)>, CsrMatrix>,
            "Newton Jacobian evaluator must return gss::CsrMatrix");

        jacobian.validate();
        if (static_cast<std::size_t>(jacobian.rows()) != state.size()) {
            reason = "Jacobian evaluator returned the wrong matrix size";
            return std::nullopt;
        }
        return jacobian;

    } catch (const std::exception& error) {
        reason = error.what();
        return std::nullopt;
    }
}

} // namespace detail

/// @brief Generic damped Newton method for R(psi) = 0.
///
/// The caller supplies two independent operations:
///   residual(state) -> R(state)
///   jacobian(state) -> J(state)
///
/// Newton solves J delta = -R with PETSc and applies
///   psi_{k+1} = psi_k + damping * delta.
///
/// SAFEGUARD, NOT GLOBALISATION. A trial state is rejected only when it cannot
/// be used: a non-finite component, a state the configured admissibility
/// predicate rejects, or a residual evaluation that fails or returns non-finite
/// values. The step is then shortened and retried. A step that merely
/// increases the residual is accepted, because requiring decrease would turn
/// this into a line search and change the method being measured.
///
/// A failed solve is reported through NewtonResult; it never terminates the
/// caller. This matters for campaigns in which a diverging state is an
/// expected outcome to be recorded, not an error to be handled.
class Newton final {
public:
    template<class ResidualEvaluator, class JacobianEvaluator, class Observer = NoNewtonObserver>
    [[nodiscard]]
    static NewtonResult solve(
        std::span<const Real> initial_state,
        ResidualEvaluator&& residual_evaluator,
        JacobianEvaluator&& jacobian_evaluator,
        PetscLinearConfig linear_config,
        NewtonSolverConfig config,
        Observer&& observer = {}) {

        using clock = std::chrono::steady_clock;
        const auto start = clock::now();

        config.validate();
        linear_config.validate();

        if (initial_state.empty()) {
            throw std::invalid_argument(
                "Newton requires a non-empty initial state");
        }
        if (!detail::newton_all_finite(initial_state)) {
            throw std::invalid_argument(
                "Newton initial state contains a non-finite value");
        }

        NewtonResult result;
        result.solution.assign(initial_state.begin(), initial_state.end());

        auto initial_attempt = detail::try_evaluate_newton_residual(
            residual_evaluator,
            std::span<const Real>(result.solution));

        if (!initial_attempt.ok) {
            result.termination = NewtonTermination::EvaluationFailure;
            result.evaluation_failure_reason =
                "initial state: " + initial_attempt.reason;
            result.total_seconds =
                std::chrono::duration<double>(clock::now() - start).count();
            return result;
        }

        auto residual = std::move(initial_attempt.values);

        const Real initial_residual = detail::newton_norm2(residual);
        const Real residual_scale =
            detail::newton_relative_denominator(initial_residual);

        result.initial_residual_norm = initial_residual;
        result.final_residual_norm = initial_residual;
        result.final_relative_residual = initial_residual / residual_scale;
        result.residual_history.push_back(initial_residual);
        result.relative_residual_history.push_back(
            result.final_relative_residual);

        std::invoke(
            observer,
            NewtonIterationInfo{
                0,
                std::span<const Real>(result.solution),
                initial_residual,
                result.final_relative_residual,
                Real{0},
                Real{0},
                Real{0},
                0,
                Real{0},
                0});

        if (initial_residual <= config.absolute_residual_tolerance) {
            result.converged = true;
            result.termination = NewtonTermination::ConvergedAbsoluteResidual;
            result.total_seconds =
                std::chrono::duration<double>(clock::now() - start).count();
            return result;
        }

        std::unique_ptr<PetscLinearSystemController> linear_solver;
        std::vector<Real> rhs(initial_state.size(), Real{0});
        std::vector<Real> previous = result.solution;
        std::vector<Real> trial(initial_state.size(), Real{0});
        std::vector<Real> step(initial_state.size(), Real{0});

        const auto state_is_admissible =
            [&config](std::span<const Real> state) {
                if (!detail::newton_all_finite(state)) {
                    return false;
                }
                if (!config.admissible) {
                    return true;
                }
                return config.admissible(state);
            };

        for (Index iteration = 0;
             iteration < config.max_iterations;
             ++iteration) {

            std::string jacobian_reason;
            auto jacobian = detail::try_evaluate_newton_jacobian(
                jacobian_evaluator,
                std::span<const Real>(result.solution),
                jacobian_reason);

            if (!jacobian.has_value()) {
                result.termination = NewtonTermination::EvaluationFailure;
                result.evaluation_failure_reason =
                    "Jacobian: " + jacobian_reason;
                break;
            }

            ++result.jacobian_evaluations;

            if (!linear_solver) {
                linear_solver = std::make_unique<PetscLinearSystemController>(
                    *jacobian,
                    linear_config);
            } else {
                linear_solver->update_matrix(*jacobian);
                ++result.jacobian_updates;
            }

            for (std::size_t i = 0; i < rhs.size(); ++i) {
                rhs[i] = -residual[i];
            }

            // The Newton correction is always solved from PETSc's zero initial
            // guess. This remains valid for KSPPREONLY+LU and avoids coupling
            // the nonlinear state to the linear correction initial guess.
            const auto linear = linear_solver->solve(rhs);
            ++result.linear_solves;
            result.total_linear_iterations += linear.iterations;
            result.max_linear_true_relative_residual = std::max(
                result.max_linear_true_relative_residual,
                static_cast<Real>(linear.true_relative_residual));

            if (!linear.converged()) {
                result.termination = NewtonTermination::LinearSolveFailure;
                result.linear_failure_reason = linear.reason_string;
                break;
            }

            previous = result.solution;

            Real length = config.damping;
            Index reductions = 0;
            bool accepted = false;
            std::vector<Real> accepted_residual;
            std::string rejection_reason;

            while (true) {

                for (std::size_t i = 0; i < trial.size(); ++i) {
                    step[i] = length * linear.solution[i];
                    trial[i] = previous[i] + step[i];
                }

                if (state_is_admissible(std::span<const Real>(trial))) {

                    auto attempt = detail::try_evaluate_newton_residual(
                        residual_evaluator,
                        std::span<const Real>(trial));

                    if (attempt.ok) {
                        accepted_residual = std::move(attempt.values);
                        accepted = true;
                        break;
                    }
                    rejection_reason = attempt.reason;

                } else {
                    rejection_reason = "trial state is not admissible";
                }

                if (reductions >= config.max_step_reductions) {
                    break;
                }

                length *= config.step_reduction_factor;
                ++reductions;
            }

            if (!accepted) {
                result.termination = NewtonTermination::StepRejected;
                result.evaluation_failure_reason = rejection_reason;
                break;
            }

            result.solution = trial;
            residual = std::move(accepted_residual);

            result.step_reductions += reductions;

            const Real step_norm = detail::newton_norm2(step);
            const Real state_norm = detail::newton_norm2(result.solution);
            const Real relative_step =
                step_norm /
                std::max(
                    Real{1},
                    detail::newton_relative_denominator(state_norm));

            const Real residual_norm = detail::newton_norm2(residual);
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
            result.damping_history.push_back(length);
            result.step_reduction_history.push_back(reductions);
            result.linear_iterations_history.push_back(linear.iterations);
            result.linear_true_relative_residual_history.push_back(
                static_cast<Real>(linear.true_relative_residual));

            std::invoke(
                observer,
                NewtonIterationInfo{
                    iteration + 1,
                    std::span<const Real>(result.solution),
                    residual_norm,
                    relative_residual,
                    step_norm,
                    relative_step,
                    length,
                    linear.iterations,
                    static_cast<Real>(linear.true_relative_residual),
                    reductions});

            if (residual_norm <= config.absolute_residual_tolerance) {
                result.converged = true;
                result.termination = NewtonTermination::ConvergedAbsoluteResidual;
                break;
            }
            if (relative_residual <= config.relative_residual_tolerance) {
                result.converged = true;
                result.termination = NewtonTermination::ConvergedRelativeResidual;
                break;
            }
        }

        result.total_seconds =
            std::chrono::duration<double>(clock::now() - start).count();
        return result;
    }
};

} // namespace gss