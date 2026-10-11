#pragma once

#include <string>
#include <vector>

#include "gsmaker/core/Types.hpp"

namespace gss {

enum class PicardTermination {
    ConvergedAbsoluteResidual,
    ConvergedRelativeResidual,
    MaximumIterations,
    LinearSolveFailure
};

[[nodiscard]]
inline const char* to_string(PicardTermination value) noexcept {
    switch (value) {
    case PicardTermination::ConvergedAbsoluteResidual:
        return "converged_absolute_residual";
    case PicardTermination::ConvergedRelativeResidual:
        return "converged_relative_residual";
    case PicardTermination::MaximumIterations:
        return "maximum_iterations";
    case PicardTermination::LinearSolveFailure:
        return "linear_solve_failure";
    }
    return "unknown";
}

/// @brief Diagnostics returned by one Picard solve.
///
/// residual_history and relative_residual_history include the initial state at
/// index zero. Step and linear histories contain one entry per completed
/// nonlinear iteration.
struct PicardResult {
    std::vector<Real> solution;

    std::vector<Real> residual_history;
    std::vector<Real> relative_residual_history;
    std::vector<Real> step_norm_history;
    std::vector<Real> relative_step_history;
    std::vector<Index> linear_iterations_history;
    std::vector<Real> linear_true_relative_residual_history;

    bool converged{false};
    PicardTermination termination{PicardTermination::MaximumIterations};
    std::string linear_failure_reason;

    Index iterations{0};
    Index linear_solves{0};
    Index total_linear_iterations{0};

    Real initial_residual_norm{0};
    Real final_residual_norm{0};
    Real final_relative_residual{0};
    Real final_step_norm{0};
    Real final_relative_step{0};
    Real max_linear_true_relative_residual{0};

    double total_seconds{0.0};
};

} // namespace gss
