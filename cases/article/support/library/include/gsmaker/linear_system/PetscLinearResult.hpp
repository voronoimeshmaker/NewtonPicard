#pragma once

#include <string>
#include <vector>

#include <petscksp.h>

#include "gsmaker/core/Types.hpp"


namespace gss {

/// @brief Result and diagnostics from one PETSc linear solve.
struct PetscLinearResult {
    std::vector<Real> solution;

    Index iterations{0};

    PetscReal ksp_residual_norm{0};
    PetscReal true_residual_norm{0};
    PetscReal rhs_norm{0};
    PetscReal true_relative_residual{0};

    KSPConvergedReason reason{
        KSP_CONVERGED_ITERATING};

    std::string reason_string;
    std::string ksp_type;
    std::string pc_type;

    double setup_seconds{0.0};
    double solve_seconds{0.0};

    [[nodiscard]]
    bool converged() const noexcept {
        return reason > 0;
    }
};

} // namespace gss
