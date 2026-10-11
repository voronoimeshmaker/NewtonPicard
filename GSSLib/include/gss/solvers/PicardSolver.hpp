#pragma once

#include <chrono>
#include <cmath>
#include <span>

#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/assembly/RHSAssembler.hpp"
#include "gss/analysis/ErrorNorms.hpp"
#include "gss/linear/PETScLinearSystem.hpp"
#include "gss/solvers/Traits.hpp"

namespace gss {

/// @brief Picard fixed-point nonlinear solver.
template<typename LinearBackend, typename RHS = RHSAssembler<>>
requires LinearSolverConcept<LinearBackend>
class PicardSolver {
public:
    /// @brief Construct from immutable configuration.
    explicit PicardSolver(const Config& config) : config_(config), backend_(config), rhs_(config) {}

    /// @brief Solve A psi = B(psi) with fixed matrix and fixed Cholesky factors.
    [[nodiscard]] SolveResult solve(MeshAdapter& mesh, std::span<const Real> initial_guess) {
        using clock = std::chrono::steady_clock;
        SolveResult result{};
        result.psi.assign(initial_guess.begin(), initial_guess.end());
        mesh.data().psi = result.psi;
        mesh.data().psi_old = result.psi;

        const auto total_start = clock::now();
        const auto ma0 = clock::now();
        Mat A = MatrixAssembler::assemble(mesh);
        // Optional fixed boundary operator supplied by a boundary-aware RHS.
        if constexpr (requires { rhs_.add_matrix_terms(mesh, A); }) {
            rhs_.add_matrix_terms(mesh, A);
        }
        if (config_.solver.diagonal_shift != Real{0}) {
            petsc_check(MatShift(A, config_.solver.diagonal_shift), "shift Picard operator");
        }
        result.timing.matrix_assembly = elapsed(ma0, clock::now());

        const auto fac0 = clock::now();
        backend_.factorize(A);
        result.timing.cholesky_factorization = elapsed(fac0, clock::now());
        result.timing.linear_setup = result.timing.cholesky_factorization;

        Vec x = PETScLinearSystem::make_vec(initial_guess);
        Vec x_tilde = PETScLinearSystem::make_vec(mesh.n_cells());
        Vec b = PETScLinearSystem::make_vec(mesh.n_cells());
        Vec residual = PETScLinearSystem::make_vec(mesh.n_cells());

        rhs_.assemble_into(mesh, x, b);
        petsc_check(MatMult(A, x, residual), "compute initial residual linear part");
        petsc_check(VecAXPY(residual, -1.0, b), "compute initial nonlinear residual");
        Real rho0 = 0.0;
        petsc_check(VecNorm(residual, NORM_2, &rho0), "compute initial nonlinear residual norm");
        rho0 = std::max(rho0, config_.solver.tolerance_absolute);

        for (Index iter = 0; iter < config_.solver.max_iter; ++iter) {
            const auto rhs0 = clock::now();
            rhs_.assemble_into(mesh, x, b);
            result.timing.rhs_assembly += elapsed(rhs0, clock::now());

            const auto solve0 = clock::now();
            backend_.solve(b, x_tilde);
            if (config_.solver.relaxation_omega != Real{1}) {
                petsc_check(VecAXPBY(x, config_.solver.relaxation_omega, Real{1} - config_.solver.relaxation_omega, x_tilde),
                            "relax Picard iterate");
            } else {
                petsc_check(VecCopy(x_tilde, x), "accept Picard iterate");
            }
            result.timing.triangular_solve += elapsed(solve0, clock::now());
            result.timing.linear_solve = result.timing.triangular_solve;

            const auto conv0 = clock::now();
            rhs_.assemble_into(mesh, x, b);
            petsc_check(MatMult(A, x, residual), "compute linear residual part");
            petsc_check(VecAXPY(residual, -1.0, b), "subtract nonlinear RHS from residual");
            Real rho = 0.0;
            petsc_check(VecNorm(residual, NORM_2, &rho), "compute nonlinear residual norm");
            result.residual_history.push_back(rho);
            result.iterations = iter + 1;
            result.timing.convergence_check += elapsed(conv0, clock::now());
            if (rho <= config_.solver.tolerance_absolute
                || rho / rho0 <= config_.solver.tolerance_relative) {
                result.converged = true;
                break;
            }
        }

        result.psi = PETScLinearSystem::to_std_vector(x);
        mesh.data().psi = result.psi;
        result.timing.total = elapsed(total_start, clock::now());
        PETScLinearSystem::destroy(residual);
        PETScLinearSystem::destroy(b);
        PETScLinearSystem::destroy(x_tilde);
        PETScLinearSystem::destroy(x);
        PETScLinearSystem::destroy(A);
        return result;
    }

private:
    static double elapsed(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    }

    Config config_;
    LinearBackend backend_;
    RHS rhs_;
};

} // namespace gss
