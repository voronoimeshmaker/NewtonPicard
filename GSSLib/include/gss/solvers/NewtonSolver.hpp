#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <span>
#include <adolc/adolc.h>
#include <adolc/sparse/sparsedrivers.h>
#include <cstdlib>
#include <limits>
#include <type_traits>

#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/assembly/RHSAssembler.hpp"
#include "gss/analysis/ErrorNorms.hpp"
#include "gss/linear/PETScLinearSystem.hpp"
#include "gss/solvers/Traits.hpp"

namespace gss {

/// @brief Sparse automatic differentiation of F(psi) using ADOL-C and ColPack.
struct ADOLCJacobian {
    template<typename RHS>
    [[nodiscard]] static Mat assemble(MeshAdapter& mesh, const RHS& rhs,
                                      std::span<const Real> psi, Real diagonal_shift = Real{0}) {
        static_assert(std::is_same_v<Real, double>, "ADOL-C requires real double PETSc scalars");
        const Index n = mesh.n_cells();
        if (n <= 0 || n > std::numeric_limits<int>::max()
            || psi.size() != static_cast<std::size_t>(n)) {
            throw std::invalid_argument("Invalid ADOL-C Jacobian dimension");
        }
        if (isTaping()) {
            throw std::logic_error("Newton Jacobian cannot interrupt an active ADOL-C recording");
        }
        // Pick an unused tag so tapes belonging to other application code survive.
        std::vector<short> tags;
        cachedTraceTags(tags);
        int candidate = 1;
        while (candidate <= std::numeric_limits<short>::max()
               && std::find(tags.begin(), tags.end(), candidate) != tags.end()) ++candidate;
        if (candidate > std::numeric_limits<short>::max()) {
            throw std::runtime_error("No free ADOL-C tape tag");
        }
        const short tag = static_cast<short>(candidate);
        struct TapeResources {
            short tag;
            bool recording{false};
            unsigned int* rows{nullptr};
            unsigned int* cols{nullptr};
            double* values{nullptr};
            ~TapeResources() {
                if (recording) trace_off();
                std::free(rows);
                std::free(cols);
                std::free(values);
                removeTape(tag, ADOLC_REMOVE_COMPLETELY);
            }
        } tape{tag};
        // Retape at every Newton state to follow state-dependent branches.
        trace_on(tag);
        tape.recording = true;
        std::vector<adouble> active(static_cast<std::size_t>(n));
        for (Index j = 0; j < n; ++j) active[j] <<= psi[j];
        const std::span<const adouble> x(active);
        for (Index i = 0; i < n; ++i) {
            adouble residual = MatrixAssembler::stencil_at(mesh, x, i)
                + diagonal_shift * x[i] - rhs.value_at(mesh, x, i);
            double value{};
            residual >>= value;
        }
        trace_off();
        tape.recording = false;
        int nnz = 0;
        int options[4] = {0, 0, 0, 0};
        const int status = sparse_jac(tag, static_cast<int>(n), static_cast<int>(n), 0,
                                      psi.data(), &nnz, &tape.rows, &tape.cols, &tape.values, options);
        if (status < 0) throw std::runtime_error("ADOL-C sparse Jacobian evaluation failed");
        Mat J{};
        petsc_check(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                n, n, 0, nullptr, 0, nullptr, &J), "create ADOL-C Jacobian");
        try {
            petsc_check(MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE),
                        "allow sparse Jacobian entries");
            Index begin{}, end{};
            petsc_check(MatGetOwnershipRange(J, &begin, &end), "get Jacobian owned rows");
            for (int k = 0; k < nnz; ++k) {
                const Index row = static_cast<Index>(tape.rows[k]);
                const Index col = static_cast<Index>(tape.cols[k]);
                if (!std::isfinite(tape.values[k])) throw std::runtime_error("Non-finite ADOL-C derivative");
                if (row >= begin && row < end) {
                    petsc_check(MatSetValue(J, row, col, tape.values[k], INSERT_VALUES),
                                "set ADOL-C Jacobian entry");
                }
            }
            petsc_check(MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY), "begin ADOL-C Jacobian assembly");
            petsc_check(MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY), "end ADOL-C Jacobian assembly");
        } catch (...) {
            MatDestroy(&J);
            throw;
        }
        return J;
    }
};

/// @brief Newton nonlinear solver with policy-selected Jacobian.
template<typename LinearBackend, typename JacobianPolicy = ADOLCJacobian, typename RHS = RHSAssembler<>>
requires LinearSolverConcept<LinearBackend>
class NewtonSolver {
public:
    /// @brief Construct from immutable configuration.
    explicit NewtonSolver(const Config& config) : config_(config), backend_(config), rhs_(config) {
        if (config.newton.jacobian != "adolc") {
            throw std::invalid_argument("newton.jacobian must be adolc; finite differences are no longer supported");
        }
    }

    /// @brief Solve F(psi) = (A + diagonal_shift I) psi - B(psi) = 0.
    [[nodiscard]] SolveResult solve(MeshAdapter& mesh, std::span<const Real> initial_guess) {
        using clock = std::chrono::steady_clock;
        SolveResult result{};
        result.psi.assign(initial_guess.begin(), initial_guess.end());
        mesh.data().psi = result.psi;
        const auto total_start = clock::now();
        Real rho0 = 0.0;
        for (Index iter = 0; iter < config_.solver.max_iter; ++iter) {
            const auto rhs0 = clock::now();
            const auto b = rhs_.assemble(mesh, result.psi);
            result.timing.rhs_assembly += elapsed(rhs0, clock::now());
            const auto Apsi = MatrixAssembler::apply_stencil(mesh, result.psi);
            std::vector<Real> F(Apsi.size(), 0.0);
            for (std::size_t i = 0; i < F.size(); ++i) {
                F[i] = Apsi[i] + config_.solver.diagonal_shift * result.psi[i] - b[i];
            }
            const Real rho = compute_norms(F, mesh.data().volume).L2;
            if (iter == 0) {
                rho0 = std::max(rho, config_.solver.tolerance_absolute);
            }
            result.residual_history.push_back(rho);
            result.iterations = iter + 1;
            if (rho <= config_.solver.tolerance_absolute || rho / rho0 <= config_.solver.tolerance_relative) {
                result.converged = true;
                break;
            }

            const auto ma0 = clock::now();
            Mat J = JacobianPolicy::assemble(mesh, rhs_, result.psi, config_.solver.diagonal_shift);
            result.timing.matrix_assembly += elapsed(ma0, clock::now());
            const auto fac0 = clock::now();
            backend_.factorize(J);
            result.timing.cholesky_factorization += elapsed(fac0, clock::now());

            for (auto& v : F) {
                v = -v;
            }
            Vec bvec = PETScLinearSystem::make_vec(F);
            std::vector<Real> delta(result.psi.size(), 0.0);
            Vec xvec = PETScLinearSystem::make_vec(delta);
            const auto solve0 = clock::now();
            backend_.solve(bvec, xvec);
            result.timing.triangular_solve += elapsed(solve0, clock::now());
            delta = PETScLinearSystem::to_std_vector(xvec);
            PETScLinearSystem::destroy(bvec);
            PETScLinearSystem::destroy(xvec);
            PETScLinearSystem::destroy(J);

            Real alpha = 1.0;
            if (config_.newton.linesearch) {
                alpha = line_search(mesh, result.psi, delta, rho);
            }
            for (std::size_t i = 0; i < result.psi.size(); ++i) {
                result.psi[i] += alpha * delta[i];
            }
            mesh.data().psi = result.psi;
        }
        result.timing.total = elapsed(total_start, clock::now());
        return result;
    }

private:
    [[nodiscard]] Real line_search(MeshAdapter& mesh, const std::vector<Real>& psi, const std::vector<Real>& delta, Real rho) const {
        Real alpha = 1.0;
        std::vector<Real> trial(psi.size());
        for (Index ls = 0; ls < config_.newton.linesearch_max_iter; ++ls) {
            for (std::size_t i = 0; i < psi.size(); ++i) {
                trial[i] = psi[i] + alpha * delta[i];
            }
            const auto b = rhs_.assemble(mesh, trial);
            const auto Apsi = MatrixAssembler::apply_stencil(mesh, trial);
            std::vector<Real> F(Apsi.size(), 0.0);
            for (std::size_t i = 0; i < F.size(); ++i) {
                F[i] = Apsi[i] + config_.solver.diagonal_shift * trial[i] - b[i];
            }
            if (compute_norms(F, mesh.data().volume).L2 <= rho) {
                return alpha;
            }
            alpha *= config_.newton.linesearch_reduction;
        }
        return alpha;
    }

    static double elapsed(std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    }

    Config config_;
    LinearBackend backend_;
    RHS rhs_;
};

} // namespace gss
