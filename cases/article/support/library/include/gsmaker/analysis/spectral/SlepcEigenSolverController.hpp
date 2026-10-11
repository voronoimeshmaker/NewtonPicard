#pragma once

#ifdef GSMAKER_USE_SLEPC

#include <slepceps.h>

#include "gsmaker/analysis/spectral/SlepcEigenConfig.hpp"
#include "gsmaker/analysis/spectral/SlepcEigenResult.hpp"
#include "gsmaker/assembly/CsrMatrix.hpp"


namespace gss {

/// @brief Reusable sequential SLEPc EPS controller for explicit CSR matrices.
///
/// This first spectral gate intentionally supports explicit matrices only.
/// Matrix-free MatShell operators are introduced separately after the EPS
/// integration itself has been verified.
class SlepcEigenSolverController final {
public:
    SlepcEigenSolverController(
        const CsrMatrix& matrix,
        SlepcEigenConfig config);

    ~SlepcEigenSolverController() noexcept;

    SlepcEigenSolverController(
        const SlepcEigenSolverController&) = delete;

    SlepcEigenSolverController& operator=(
        const SlepcEigenSolverController&) = delete;

    SlepcEigenSolverController(
        SlepcEigenSolverController&&) = delete;

    SlepcEigenSolverController& operator=(
        SlepcEigenSolverController&&) = delete;

    [[nodiscard]]
    SlepcEigenResult solve();

private:
    void configure_eps();

    Mat matrix_{nullptr};
    EPS eps_{nullptr};

    SlepcEigenConfig config_{};
};

} // namespace gss

#endif // GSMAKER_USE_SLEPC
