#pragma once

#include <array>
#include <string>
#include <vector>

#include "gss/core/Types.hpp"

namespace gss {

/// @brief Mesh generation and domain parameters.
struct MeshConfig {
    std::string type{"cartesian_square"};
    Index n_generators{1024};
    Index n_per_direction{16};
    Index n_per_direction_y{0}; // 0 -> use n_per_direction for both
    Index random_seed{12345};
    std::array<Real, 4> domain{0.0, 1.0, 0.0, 1.0};
    bool hilbert_reorder{true};
    /// @brief Cell renumbering applied after mesh construction.
    /// natural (identity) | random | morton | hilbert.
    std::string reordering{"natural"};
};

/// @brief Nonlinear and linear solver runtime parameters.
struct SolverConfig {
    std::string nonlinear{"picard"};
    std::string linear{"cholesky_sparse"};
    Index max_iter{100};
    Real tolerance_relative{1.0e-8};
    Real tolerance_absolute{1.0e-12};
    Real relaxation_omega{1.0};
    Real rhs_constant{1.0};
    /// @brief Fixed denominator regularisation for production dF/dPsi recovery.
    Real source_regularisation_delta{1.0e-10};
    Real diagonal_shift{1.0};
    std::string initial_guess{"zero"};
    Real perturbation_amplitude{0.0};
    std::string perturbation_type{"gaussian"};
};

/// @brief Newton method parameters.
struct NewtonConfig {
    std::string jacobian{"adolc"};
    bool linesearch{true};
    Index linesearch_max_iter{10};
    Real linesearch_reduction{0.5};
};

/// @brief PETSc KSP and PC parameters.
struct PETScConfig {
    std::string ksp_type{"preonly"};
    std::string pc_type{"cholesky"};
    std::string mat_solver_type{"mumps"};
    std::string mumps_reorder{"pord"};
    Index mumps_icntl_14{20};
};

/// @brief Analytical benchmark parameters.
struct BenchmarkConfig {
    std::string type{"fadeev"};
    Real fadeev_epsilon{0.1};
    Real fadeev_psi0{1.0};
    Real fadeev_k{1.0};
    /// @brief Shape parameter for the rotated VoroGradToolkit Fadeev benchmark.
    Real fadeev_alpha{0.225};
    /// @brief Rigid rotation of the Fadeev benchmark, in degrees.
    Real fadeev_theta_degrees{5.7};
    bool compute_error{true};
    std::vector<std::string> error_norms{"L1", "L2", "Linf"};
    Real kappa_island_threshold{100.0};
};

/// @brief Convergence study parameters.
struct ConvergenceStudyConfig {
    bool enabled{false};
    std::vector<Index> n_generators_sequence{64, 256, 1024, 4096, 16384};
    std::string output_csv{"convergence_study.csv"};
};

/// @brief Gradient extrema reporting parameters.
struct GradientAnalysisConfig {
    bool enabled{true};
    Index top_k{20};
    std::vector<std::string> quantities{"grad_psi", "grad_psi_inv", "source_term"};
    Real min_distance{2.0};
    Real grad_critical_threshold{0.01};
    std::string output_format{"both"};
    std::string filename_extrema{"gradient_extrema.csv"};
    std::string filename_critical{"critical_points.csv"};
};

/// @brief Output and timing report parameters.
struct OutputConfig {
    std::string format{"vtk"};
    bool write_solution{true};
    bool write_error{true};
    bool timing_report{true};
    bool timing_components{true};
};

/// @brief Complete application configuration loaded from YAML.
struct Config {
    MeshConfig mesh;
    SolverConfig solver;
    NewtonConfig newton;
    PETScConfig petsc;
    BenchmarkConfig benchmark;
    ConvergenceStudyConfig convergence_study;
    GradientAnalysisConfig gradient_analysis;
    OutputConfig output;
};

} // namespace gss
