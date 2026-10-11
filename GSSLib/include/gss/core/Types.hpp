#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include <petsc.h>

namespace gss {

using Real = PetscScalar;
using Index = PetscInt;

/// @brief Hot cell-centered solver data stored as structure-of-arrays.
struct alignas(64) CellData {
    std::vector<Real> psi;
    std::vector<Real> psi_old;
    std::vector<Real> source;
    std::vector<Real> volume;
    std::vector<Real> grad_psi_x;
    std::vector<Real> grad_psi_y;
    std::vector<Real> grad_F_x;
    std::vector<Real> grad_F_y;
    std::vector<Real> kappa;
    Index n_cells{0};

    /// @brief Resize all hot arrays to the requested number of cells.
    void resize(Index n);
};

/// @brief Cold Voronoi topology in CSR format.
struct MeshTopology {
    std::vector<Index> neighbor_start;
    std::vector<Index> neighbor_list;
    std::vector<Real> interface_area;
    std::vector<Real> distance;
    std::vector<Real> coord_x;
    std::vector<Real> coord_y;
    std::vector<unsigned char> is_boundary_cell;
};

/// @brief Wall-clock timing breakdown for a nonlinear solve.
struct TimingReport {
    double matrix_assembly{0};
    double cholesky_factorization{0};
    double linear_setup{0};
    double gradient_reconstruction{0};
    double rhs_assembly{0};
    double triangular_solve{0};
    double linear_solve{0};
    double convergence_check{0};
    double error_computation{0};
    double gradient_extrema{0};
    double total{0};
};

/// @brief Nonlinear solve result including solution and residual history.
struct SolveResult {
    std::vector<Real> psi;
    std::vector<Real> residual_history;
    TimingReport timing;
    bool converged{false};
    Index iterations{0};
};

} // namespace gss
