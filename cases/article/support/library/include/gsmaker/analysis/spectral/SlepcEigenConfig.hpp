#pragma once

#ifdef GSMAKER_USE_SLEPC

#include <cmath>
#include <stdexcept>
#include <string>

#include <slepceps.h>


namespace gss {

/// @brief Runtime configuration for one SLEPc EPS eigensolver.
///
/// EPS solver names remain native SLEPc strings.  GSMakerLib does not mirror
/// the SLEPc solver catalogue with a parallel enum/switch hierarchy.
struct SlepcEigenConfig {
    std::string eps_type{"krylovschur"};
    std::string options_prefix{};

    PetscInt requested_eigenpairs{2};
    PetscReal tolerance{1.0e-10};
    PetscInt max_iterations{1000};

    EPSProblemType problem_type{EPS_NHEP};
    EPSWhich which{EPS_LARGEST_MAGNITUDE};

    bool set_from_options{true};
    bool compute_eigenvectors{true};

    void validate() const {
        if (eps_type.empty()) {
            throw std::invalid_argument(
                "SlepcEigenConfig eps_type must not be empty");
        }

        if (requested_eigenpairs <= 0) {
            throw std::invalid_argument(
                "SlepcEigenConfig requested_eigenpairs must be positive");
        }

        if (!std::isfinite(static_cast<double>(tolerance)) ||
            tolerance <= PetscReal{0}) {
            throw std::invalid_argument(
                "SlepcEigenConfig tolerance must be finite and positive");
        }

        if (max_iterations <= 0) {
            throw std::invalid_argument(
                "SlepcEigenConfig max_iterations must be positive");
        }
    }
};

} // namespace gss

#endif // GSMAKER_USE_SLEPC
