#pragma once

#include <stdexcept>

#include <petscksp.h>

#include "gss/core/Config.hpp"
#include "gss/linear/PETScLinearSystem.hpp"

namespace gss {

/// @brief PETSc KSP backend with reusable operators and preconditioner setup.
class CholeskyBackend {
public:
    /// @brief Construct with immutable PETSc settings.
    explicit CholeskyBackend(const Config& config) : config_(config) {}

    /// @brief Destroy owned KSP object.
    ~CholeskyBackend();

    CholeskyBackend(const CholeskyBackend&) = delete;
    CholeskyBackend& operator=(const CholeskyBackend&) = delete;
    CholeskyBackend(CholeskyBackend&& other) noexcept;
    CholeskyBackend& operator=(CholeskyBackend&& other) noexcept;

    /// @brief Configure and set up the KSP/PC for the provided matrix.
    void factorize(Mat A);

    /// @brief Solve the factorized system for x.
    void solve(Vec b, Vec x);

    /// @brief Last KSPSolve iteration count.
    [[nodiscard]] Index last_iterations() const noexcept { return last_iterations_; }

    /// @brief Borrow the configured KSP for read-only diagnostic queries.
    [[nodiscard]] KSP native_ksp() const noexcept { return ksp_; }

private:
    Config config_;
    KSP ksp_{nullptr};
    Index last_iterations_{0};
};

using PETScKSPBackend = CholeskyBackend;

} // namespace gss
