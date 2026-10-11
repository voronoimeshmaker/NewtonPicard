#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include <petscksp.h>

#include "gsmaker/core/Config.hpp"


namespace gss {

/// @brief Runtime configuration for one reusable PETSc KSP controller.
///
/// Solver/preconditioner names are PETSc strings.  GSMakerLib deliberately
/// does not maintain a parallel enum or switch table for KSP/PC types.
struct PetscLinearConfig {
    std::string ksp_type{"cg"};
    std::string pc_type{"gamg"};
    std::string mat_solver_type{};
    std::string options_prefix{};

    PetscReal relative_tolerance{1.0e-12};
    PetscReal absolute_tolerance{1.0e-14};
    PetscReal divergence_tolerance{1.0e8};

    PetscInt max_iterations{5000};

    bool set_from_options{true};

    /// Mathematical hints supplied by the caller; they are not inferred here.
    bool matrix_symmetric_hint{false};
    bool matrix_spd_hint{false};

    [[nodiscard]]
    static PetscLinearConfig from_config(
        const Config& config) {

        PetscLinearConfig result;

        result.ksp_type =
            config.petsc.ksp_type;

        result.pc_type =
            config.petsc.pc_type;

        result.mat_solver_type =
            config.petsc.mat_solver_type;

        result.relative_tolerance =
            static_cast<PetscReal>(
                config.solver.tolerance_relative);

        result.absolute_tolerance =
            static_cast<PetscReal>(
                config.solver.tolerance_absolute);

        result.max_iterations =
            static_cast<PetscInt>(
                config.solver.max_iter);

        result.validate();

        return result;
    }

    void validate() const {

        if (ksp_type.empty()) {
            throw std::invalid_argument(
                "PetscLinearConfig ksp_type must not be empty");
        }

        if (pc_type.empty()) {
            throw std::invalid_argument(
                "PetscLinearConfig pc_type must not be empty");
        }

        if (!std::isfinite(
                static_cast<double>(
                    relative_tolerance)) ||
            relative_tolerance < PetscReal{0}) {

            throw std::invalid_argument(
                "PetscLinearConfig relative tolerance must be finite and non-negative");
        }

        if (!std::isfinite(
                static_cast<double>(
                    absolute_tolerance)) ||
            absolute_tolerance < PetscReal{0}) {

            throw std::invalid_argument(
                "PetscLinearConfig absolute tolerance must be finite and non-negative");
        }

        if (!std::isfinite(
                static_cast<double>(
                    divergence_tolerance)) ||
            divergence_tolerance <= PetscReal{1}) {

            throw std::invalid_argument(
                "PetscLinearConfig divergence tolerance must be finite and greater than one");
        }

        if (max_iterations <= 0) {
            throw std::invalid_argument(
                "PetscLinearConfig max_iterations must be positive");
        }
    }
};

} // namespace gss
