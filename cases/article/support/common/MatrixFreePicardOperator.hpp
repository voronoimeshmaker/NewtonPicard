#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include <petscksp.h>
#include <slepceps.h>

#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/core/Types.hpp"

namespace study {
namespace spectral {

/// @brief Leading eigenvalues of the Picard operator without forming it.
///
/// WHY THIS EXISTS.  study::build_explicit_G forms G = A^{-1} B one column at a
/// time: n linear solves and n^2 stored entries, followed by a dense eigenvalue
/// problem.  That is the second O(n^3) cost in the campaigns, independent of
/// the Jacobian assembly, and it is the one that binds the spectral studies,
/// which are the core of the manuscript.
///
/// The manuscript already states the remedy: the action of G needs only an
/// application of B followed by a solve with A, so a Krylov eigensolver can
/// work on the operator directly.  A is invariant across the whole many-query
/// calculation, so its factorisation is computed once here and reused by every
/// matrix-vector product; that reuse is the amortisation argument of the paper,
/// exercised rather than asserted.
///
/// THREE OPERATORS, ONE SHELL.  The campaigns currently form three different
/// dense objects.  All three are the same action with a different final
/// combination, so they are modes here:
///
///   Picard                 y = A^{-1} B x                  spectrum gives rho
///   RelaxedPicard          y = (1 - w) x + w A^{-1} B x     stability under
///                                                           relaxation
///   PreconditionedJacobian y = x - A^{-1} B x               spectrum of
///                                                           A^{-1} J
///
/// COST IS REPORTED, NOT ASSUMED.  The number of operator applications and the
/// factorisation time are returned separately from the eigensolve time, because
/// the manuscript promises to measure what computing rho costs.
enum class OperatorMode {
    Picard,
    RelaxedPicard,
    PreconditionedJacobian
};

namespace detail {

struct ShellContext {
    Mat b_matrix{nullptr};
    KSP a_solver{nullptr};
    Vec work{nullptr};
    PetscScalar omega{1.0};
    OperatorMode mode{OperatorMode::Picard};
    PetscInt applications{0};
};

inline PetscErrorCode shell_multiply(Mat shell, Vec x, Vec y) {
    PetscFunctionBeginUser;

    ShellContext* context = nullptr;
    PetscCall(MatShellGetContext(shell, &context));

    ++context->applications;

    PetscCall(MatMult(context->b_matrix, x, context->work));
    PetscCall(KSPSolve(context->a_solver, context->work, y));

    if (context->mode == OperatorMode::RelaxedPicard) {
        // y <- (1 - w) x + w y
        const PetscScalar one{1.0};
        PetscCall(VecAXPBY(y, one - context->omega, context->omega, x));
    } else if (context->mode == OperatorMode::PreconditionedJacobian) {
        // y <- x - y
        const PetscScalar one{1.0};
        const PetscScalar minus_one{-1.0};
        PetscCall(VecAXPBY(y, one, minus_one, x));
    }

    PetscFunctionReturn(PETSC_SUCCESS);
}

inline void check(PetscErrorCode code, const char* what) {
    if (code != PETSC_SUCCESS) {
        throw std::runtime_error(
            std::string{"MatrixFreePicardOperator: "} + what + " failed");
    }
}

} // namespace detail

/// Owns the factorisation of A, a PETSc copy of B, and the shell operator.
///
/// The CSR arrays of B are copied into members because
/// MatCreateSeqAIJWithArrays does not take ownership and requires them to
/// outlive the matrix.
class MatrixFreePicardOperator final {
public:
    MatrixFreePicardOperator(
        const gss::CsrMatrix& a,
        const gss::CsrMatrix& b,
        OperatorMode mode = OperatorMode::Picard,
        gss::Real omega = gss::Real{1},
        bool a_is_spd = true) {

        const auto rows = static_cast<PetscInt>(a.rows());

        if (rows <= 0) {
            throw std::invalid_argument("MatrixFreePicardOperator requires a non-empty matrix");
        }
        if (b.rows() != a.rows()) {
            throw std::invalid_argument("MatrixFreePicardOperator A and B sizes differ");
        }

        cells_ = static_cast<std::size_t>(rows);

        const auto a_row = a.row_start();
        const auto a_col = a.column();
        const auto a_val = a.value();

        a_row_.assign(a_row.begin(), a_row.end());
        a_col_.assign(a_col.begin(), a_col.end());
        a_val_.assign(a_val.begin(), a_val.end());

        const auto b_row = b.row_start();
        const auto b_col = b.column();
        const auto b_val = b.value();

        b_row_.assign(b_row.begin(), b_row.end());
        b_col_.assign(b_col.begin(), b_col.end());
        b_val_.assign(b_val.begin(), b_val.end());

        detail::check(
            MatCreateSeqAIJWithArrays(
                PETSC_COMM_SELF, rows, rows,
                a_row_.data(), a_col_.data(), a_val_.data(), &a_matrix_),
            "creating A");

        detail::check(
            MatCreateSeqAIJWithArrays(
                PETSC_COMM_SELF, rows, rows,
                b_row_.data(), b_col_.data(), b_val_.data(), &context_.b_matrix),
            "creating B");

        // A is invariant over the whole many-query calculation, so it is
        // factorised once here and reused by every operator application.
        const auto factorisation_start = std::chrono::steady_clock::now();

        detail::check(KSPCreate(PETSC_COMM_SELF, &context_.a_solver), "creating the A solver");
        detail::check(KSPSetOperators(context_.a_solver, a_matrix_, a_matrix_), "setting A");
        detail::check(KSPSetType(context_.a_solver, KSPPREONLY), "selecting preonly");

        PC preconditioner = nullptr;
        detail::check(KSPGetPC(context_.a_solver, &preconditioner), "retrieving the preconditioner");
        detail::check(
            PCSetType(preconditioner, a_is_spd ? PCCHOLESKY : PCLU),
            "selecting the factorisation");

        detail::check(KSPSetUp(context_.a_solver), "factorising A");

        factorisation_seconds_ =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - factorisation_start).count();

        detail::check(
            MatCreateVecs(context_.b_matrix, nullptr, &context_.work),
            "creating the work vector");

        context_.omega = static_cast<PetscScalar>(omega);
        context_.mode = mode;
        context_.applications = 0;

        detail::check(
            MatCreateShell(
                PETSC_COMM_SELF, rows, rows, rows, rows, &context_, &shell_),
            "creating the shell operator");

        detail::check(
            MatShellSetOperation(
                shell_, MATOP_MULT,
                reinterpret_cast<void (*)(void)>(detail::shell_multiply)),
            "installing the shell product");
    }

    MatrixFreePicardOperator(const MatrixFreePicardOperator&) = delete;
    MatrixFreePicardOperator& operator=(const MatrixFreePicardOperator&) = delete;
    MatrixFreePicardOperator(MatrixFreePicardOperator&&) = delete;
    MatrixFreePicardOperator& operator=(MatrixFreePicardOperator&&) = delete;

    ~MatrixFreePicardOperator() {
        // Destructors must not throw, so return codes are deliberately ignored.
        if (shell_ != nullptr) { MatDestroy(&shell_); }
        if (context_.work != nullptr) { VecDestroy(&context_.work); }
        if (context_.a_solver != nullptr) { KSPDestroy(&context_.a_solver); }
        if (context_.b_matrix != nullptr) { MatDestroy(&context_.b_matrix); }
        if (a_matrix_ != nullptr) { MatDestroy(&a_matrix_); }
    }

    /// Replace the values of B without touching the factorisation of A.
    ///
    /// WHY THIS EXISTS.  A is invariant across the whole many-query
    /// calculation and B is the only operator that changes with the free
    /// function.  Constructing a new operator per case would refactorise A
    /// every time, which is the one cost the manuscript argues is paid once.
    /// The sparsity pattern of B is the face-neighbour graph and therefore
    /// fixed by the mesh, so a new case changes the values alone.
    ///
    /// MatCreateSeqAIJWithArrays shares b_val_ with PETSc rather than copying
    /// it, so the update is a copy into the member followed by a reassembly.
    void update_b_values(std::span<const gss::Real> values) {
        if (values.size() != b_val_.size()) {
            throw std::invalid_argument(
                "MatrixFreePicardOperator: the number of entries in B changed; "
                "its pattern is mesh-invariant and must not");
        }
        std::copy(values.begin(), values.end(), b_val_.begin());
        detail::check(
            MatAssemblyBegin(context_.b_matrix, MAT_FINAL_ASSEMBLY),
            "reassembling B after a value update");
        detail::check(
            MatAssemblyEnd(context_.b_matrix, MAT_FINAL_ASSEMBLY),
            "reassembling B after a value update");
    }

    /// Change the relaxation factor used by OperatorMode::RelaxedPicard.
    /// Sweeping omega needs no new factorisation and no new operator.
    void set_omega(gss::Real omega) noexcept {
        context_.omega = static_cast<PetscScalar>(omega);
    }

    /// Switch between the three spectra the campaigns need.  The factorisation
    /// of A and the copy of B are shared by all three.
    void set_mode(OperatorMode mode) noexcept { context_.mode = mode; }

    [[nodiscard]] OperatorMode mode() const noexcept { return context_.mode; }
    [[nodiscard]] gss::Real omega() const noexcept {
        return static_cast<gss::Real>(PetscRealPart(context_.omega));
    }

    [[nodiscard]] Mat mat() const noexcept { return shell_; }
    [[nodiscard]] std::size_t cells() const noexcept { return cells_; }
    [[nodiscard]] std::size_t b_entries() const noexcept { return b_val_.size(); }
    [[nodiscard]] PetscInt applications() const noexcept { return context_.applications; }
    [[nodiscard]] double factorisation_seconds() const noexcept { return factorisation_seconds_; }

    void reset_application_count() noexcept { context_.applications = 0; }

private:
    std::size_t cells_{};

    std::vector<PetscInt> a_row_;
    std::vector<PetscInt> a_col_;
    std::vector<PetscScalar> a_val_;

    std::vector<PetscInt> b_row_;
    std::vector<PetscInt> b_col_;
    std::vector<PetscScalar> b_val_;

    Mat a_matrix_{nullptr};
    Mat shell_{nullptr};

    detail::ShellContext context_{};

    double factorisation_seconds_{0.0};
};

struct MatrixFreeSpectrum {
    gss::Real rho{0};
    std::complex<gss::Real> mu1{};
    std::complex<gss::Real> mu2{};

    PetscInt converged_eigenpairs{0};
    PetscInt eps_iterations{0};
    EPSConvergedReason reason{EPS_CONVERGED_ITERATING};
    gss::Real max_relative_error{0};

    /// Operator applications used by the eigensolve.  Each one costs a sparse
    /// product with B and one triangular solve with the factorised A; this is
    /// the recurring cost of learning rho, and the figure the many-query
    /// argument needs.
    PetscInt applications{0};

    /// Time to factorise A, paid once for the whole calculation.
    double factorisation_seconds{0.0};

    /// Time inside the eigensolve, excluding the factorisation.
    double eigensolve_seconds{0.0};
};

/// Leading eigenvalues of the operator, largest magnitude first.
[[nodiscard]] inline MatrixFreeSpectrum leading_spectrum(
    MatrixFreePicardOperator& op,
    PetscInt requested_eigenpairs = 6,
    PetscReal tolerance = 1.0e-11,
    PetscInt max_iterations = 2000) {

    const auto n = static_cast<PetscInt>(op.cells());

    EPS eps = nullptr;
    detail::check(EPSCreate(PETSC_COMM_SELF, &eps), "creating the eigensolver");

    MatrixFreeSpectrum summary;

    try {
        detail::check(EPSSetOperators(eps, op.mat(), nullptr), "setting the operator");
        detail::check(EPSSetProblemType(eps, EPS_NHEP), "declaring a non-Hermitian problem");
        detail::check(EPSSetWhichEigenpairs(eps, EPS_LARGEST_MAGNITUDE), "selecting the target");
        detail::check(EPSSetType(eps, EPSKRYLOVSCHUR), "selecting Krylov-Schur");

        const PetscInt nev = std::min<PetscInt>(std::max<PetscInt>(requested_eigenpairs, 1), n);

        detail::check(
            EPSSetDimensions(eps, nev, PETSC_DEFAULT, PETSC_DEFAULT),
            "setting the subspace size");
        detail::check(
            EPSSetTolerances(eps, tolerance, max_iterations),
            "setting the tolerances");

        op.reset_application_count();

        const auto solve_start = std::chrono::steady_clock::now();
        detail::check(EPSSolve(eps), "solving the eigenvalue problem");
        summary.eigensolve_seconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - solve_start).count();

        summary.applications = op.applications();
        summary.factorisation_seconds = op.factorisation_seconds();

        detail::check(EPSGetConverged(eps, &summary.converged_eigenpairs), "counting eigenpairs");
        detail::check(EPSGetIterationNumber(eps, &summary.eps_iterations), "reading the iteration count");
        detail::check(EPSGetConvergedReason(eps, &summary.reason), "reading the reason");

        if (summary.converged_eigenpairs < 1) {
            throw std::runtime_error(
                "MatrixFreePicardOperator: the eigensolver converged no eigenpair");
        }

        std::vector<std::complex<double>> eigenvalue;
        eigenvalue.reserve(static_cast<std::size_t>(summary.converged_eigenpairs));

        for (PetscInt i = 0; i < summary.converged_eigenpairs; ++i) {

            PetscScalar real_part{};
            PetscScalar imaginary_part{};

            detail::check(
                EPSGetEigenvalue(eps, i, &real_part, &imaginary_part),
                "reading an eigenvalue");

#if defined(PETSC_USE_COMPLEX)
            eigenvalue.emplace_back(
                static_cast<double>(PetscRealPart(real_part)),
                static_cast<double>(PetscImaginaryPart(real_part)));
#else
            eigenvalue.emplace_back(
                static_cast<double>(real_part),
                static_cast<double>(imaginary_part));
#endif

            PetscReal error{};
            detail::check(
                EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error),
                "computing an eigenpair error");

            summary.max_relative_error = std::max(
                summary.max_relative_error,
                static_cast<gss::Real>(error));
        }

        std::sort(
            eigenvalue.begin(),
            eigenvalue.end(),
            [](const std::complex<double>& lhs, const std::complex<double>& rhs) {
                return std::abs(lhs) > std::abs(rhs);
            });

        summary.rho = static_cast<gss::Real>(std::abs(eigenvalue.front()));
        summary.mu1 = {
            static_cast<gss::Real>(eigenvalue.front().real()),
            static_cast<gss::Real>(eigenvalue.front().imag())};

        if (eigenvalue.size() > 1U) {
            summary.mu2 = {
                static_cast<gss::Real>(eigenvalue[1].real()),
                static_cast<gss::Real>(eigenvalue[1].imag())};
        }

    } catch (...) {
        EPSDestroy(&eps);
        throw;
    }

    EPSDestroy(&eps);
    return summary;
}

} // namespace spectral
} // namespace study