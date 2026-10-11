#pragma once

#include <span>


#include <petscksp.h>

#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/core/Types.hpp"
#include "gsmaker/linear_system/PetscLinearConfig.hpp"
#include "gsmaker/linear_system/PetscLinearResult.hpp"


namespace gss {

/// @brief Reusable sequential PETSc linear-system controller.
///
/// Matrix/KSP/preconditioner setup is performed once.  Subsequent solve()
/// calls reuse the same PETSc objects, which is the required ownership model
/// for later Picard and many-query workloads with a constant elliptic operator.
class PetscLinearSystemController final {
public:
    PetscLinearSystemController(
        const CsrMatrix& matrix,
        PetscLinearConfig config);

    ~PetscLinearSystemController() noexcept;

    PetscLinearSystemController(
        const PetscLinearSystemController&) = delete;

    PetscLinearSystemController& operator=(
        const PetscLinearSystemController&) = delete;

    PetscLinearSystemController(
        PetscLinearSystemController&&) = delete;

    PetscLinearSystemController& operator=(
        PetscLinearSystemController&&) = delete;

    /// @brief Solve A x = rhs.
    ///
    /// If initial_guess is empty, PETSc starts from zero.  Otherwise it must
    /// have exactly the matrix dimension and KSP is told that the initial guess
    /// is non-zero.
    [[nodiscard]]
    PetscLinearResult solve(
        std::span<const Real> rhs,
        std::span<const Real> initial_guess = {});

    /// @brief Replace numerical matrix values while retaining controller/KSP.
    ///
    /// The matrix dimensions and sparsity pattern must be unchanged.
    void update_matrix(
        const CsrMatrix& matrix);

    [[nodiscard]]
    Index rows() const noexcept {
        return rows_;
    }

    [[nodiscard]]
    double setup_seconds() const noexcept {
        return setup_seconds_;
    }

private:
    void configure_ksp();

    Mat matrix_{nullptr};
    KSP ksp_{nullptr};

    PetscLinearConfig config_{};

    Index rows_{0};

    double setup_seconds_{0.0};
};

} // namespace gss
