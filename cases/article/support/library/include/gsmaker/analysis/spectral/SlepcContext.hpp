#pragma once

#ifdef GSMAKER_USE_SLEPC

#include <slepcsys.h>


namespace gss {

/// @brief Process-level SLEPc lifetime guard.
///
/// SlepcContext initializes SLEPc only when it has not already been
/// initialized.  SLEPc initializes PETSc/MPI when required, so a live
/// SlepcContext also satisfies GSMakerLib components that require PETSc.
/// The guard finalizes SLEPc only when it owns the SLEPc initialization.
///
/// Keep SlepcContext alive longer than every SLEPc object created by
/// GSMakerLib.
class SlepcContext final {
public:
    SlepcContext();
    ~SlepcContext() noexcept;

    SlepcContext(
        const SlepcContext&) = delete;

    SlepcContext& operator=(
        const SlepcContext&) = delete;

    SlepcContext(
        SlepcContext&&) = delete;

    SlepcContext& operator=(
        SlepcContext&&) = delete;

    [[nodiscard]]
    bool owns_initialization() const noexcept {
        return owns_initialization_;
    }

    [[nodiscard]]
    static bool initialized();

private:
    bool owns_initialization_{false};
};

} // namespace gss

#endif // GSMAKER_USE_SLEPC
