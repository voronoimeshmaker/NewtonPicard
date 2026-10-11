#pragma once

#include <string>
#include <vector>

#include "gsmaker/core/Types.hpp"

namespace gss {

enum class NewtonTermination {
    ConvergedAbsoluteResidual,
    ConvergedRelativeResidual,
    MaximumIterations,
    LinearSolveFailure,

    /// The correction could not be reduced to a step whose state was finite,
    /// admissible and evaluable within the configured number of reductions.
    StepRejected,

    /// The residual or the Jacobian could not be evaluated at a state the
    /// method had already accepted, including the initial state.
    EvaluationFailure
};

[[nodiscard]]
inline const char* to_string(NewtonTermination value) noexcept {
    switch (value) {
    case NewtonTermination::ConvergedAbsoluteResidual:
        return "converged_absolute_residual";
    case NewtonTermination::ConvergedRelativeResidual:
        return "converged_relative_residual";
    case NewtonTermination::MaximumIterations:
        return "maximum_iterations";
    case NewtonTermination::LinearSolveFailure:
        return "linear_solve_failure";
    case NewtonTermination::StepRejected:
        return "step_rejected";
    case NewtonTermination::EvaluationFailure:
        return "evaluation_failure";
    }
    return "unknown";
}

/// @brief Diagnostics from one Newton solve.
///
/// residual histories contain the initial state at index zero. Step, damping
/// and linear histories contain one entry per completed Newton correction.
struct NewtonResult {
    std::vector<Real> solution;

    std::vector<Real> residual_history;
    std::vector<Real> relative_residual_history;
    std::vector<Real> step_norm_history;
    std::vector<Real> relative_step_history;
    std::vector<Real> damping_history;
    std::vector<Index> linear_iterations_history;
    std::vector<Real> linear_true_relative_residual_history;

    /// Step reductions applied per accepted correction. A non-zero entry means
    /// the full damped step was not evaluable and had to be shortened.
    std::vector<Index> step_reduction_history;

    bool converged{false};
    NewtonTermination termination{NewtonTermination::MaximumIterations};
    std::string linear_failure_reason;

    /// Why a residual or Jacobian evaluation failed, when one did.
    std::string evaluation_failure_reason;

    Index iterations{0};
    Index jacobian_evaluations{0};
    Index jacobian_updates{0};
    Index linear_solves{0};
    Index total_linear_iterations{0};

    /// Total step reductions over the whole solve. Zero means every damped
    /// step was taken as configured.
    Index step_reductions{0};

    Real initial_residual_norm{0};
    Real final_residual_norm{0};
    Real final_relative_residual{0};
    Real final_step_norm{0};
    Real final_relative_step{0};
    Real max_linear_true_relative_residual{0};

    double total_seconds{0.0};
};

} // namespace gss