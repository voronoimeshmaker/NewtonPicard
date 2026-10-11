#pragma once

#include <array>
#include <string>

#include "gsmaker/core/Types.hpp"

namespace gss {

/// @brief Temporary mesh configuration used while mesh construction
/// is being separated from solver-facing mesh representation.
struct MeshConfig {
    std::string type{"cartesian_square"};

    Index n_generators{1024};
    Index n_per_direction{16};
    Index n_per_direction_y{0};
    Index random_seed{12345};

    std::array<Real, 4> domain{
        Real{0.0},
        Real{1.0},
        Real{0.0},
        Real{1.0}
    };

    // Transitional compatibility fields. GSMakerLib no longer performs mesh
    // reordering; non-natural ordering must be applied through
    // VoronoiMeshMaker before construction of the solver-facing mesh.
    bool hilbert_reorder{true};
    std::string reordering{"natural"};
};


/// @brief Configuration common to nonlinear solution procedures.
struct SolverConfig {
    std::string nonlinear{"picard"};
    std::string linear{"cholesky_sparse"};

    Index max_iter{100};

    Real tolerance_relative{1.0e-8};
    Real tolerance_absolute{1.0e-12};

    Real relaxation_omega{1.0};

    // Transitional fields. Their final ownership will be decided
    // when the nonlinear/source interfaces are migrated.
    Real rhs_constant{1.0};
    Real source_regularisation_delta{1.0e-10};
    Real diagonal_shift{1.0};

    std::string initial_guess{"zero"};
    Real perturbation_amplitude{0.0};
    std::string perturbation_type{"gaussian"};
};



/// @brief Newton-specific numerical configuration.
struct NewtonConfig {
    std::string jacobian{"finite_difference"};

    Real fd_epsilon{1.0e-7};

    bool linesearch{true};
    Index linesearch_max_iter{10};
    Real linesearch_reduction{0.5};
};


/// @brief PETSc linear-system configuration.
struct PETScConfig {
    std::string ksp_type{"preonly"};
    std::string pc_type{"cholesky"};
    std::string mat_solver_type{"mumps"};
    std::string mumps_reorder{"pord"};

    Index mumps_icntl_14{20};
};


/// @brief Configuration owned by GSMakerLib.
///
/// Problem-, benchmark-, study- and article-specific parameters do not
/// belong here. They are supplied through ConfigTraits specialisations.
struct CommonConfig {
    MeshConfig mesh{};
    SolverConfig solver{};
    NewtonConfig newton{};
    PETScConfig petsc{};
};


/// @brief Transitional compatibility name.
///
/// Existing GSMakerLib code may continue using gss::Config while the
/// architecture is migrated incrementally.
using Config = CommonConfig;

} // namespace gss
