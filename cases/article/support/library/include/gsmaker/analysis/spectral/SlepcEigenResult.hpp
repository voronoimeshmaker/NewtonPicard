#pragma once

#ifdef GSMAKER_USE_SLEPC

#include <complex>
#include <string>
#include <vector>


#include <slepceps.h>

namespace gss {

/// @brief One converged right eigenpair and SLEPc residual diagnostics.
struct SlepcEigenpair {
    std::complex<PetscReal> eigenvalue{};
    std::vector<std::complex<PetscReal>> eigenvector{};

    PetscReal absolute_error{0};
    PetscReal relative_error{0};
    PetscReal backward_error{0};
};


/// @brief Result of one EPS solve.
struct SlepcEigenResult {
    std::vector<SlepcEigenpair> eigenpairs{};

    PetscInt requested_eigenpairs{0};
    PetscInt converged_eigenpairs{0};
    PetscInt iterations{0};

    EPSConvergedReason reason{EPS_CONVERGED_ITERATING};
    std::string eps_type{};

    double solve_seconds{0.0};

    [[nodiscard]]
    bool converged_requested() const noexcept {
        return reason > 0 &&
            converged_eigenpairs >= requested_eigenpairs;
    }
};

} // namespace gss

#endif // GSMAKER_USE_SLEPC
