#pragma once

#include <string_view>

#include <petscsys.h>

namespace gss {

/// @brief Convert a non-zero PETSc error code into a C++ exception.
void petsc_check(
    PetscErrorCode error,
    std::string_view action);


/// @brief Process-level PETSc lifetime guard.
///
/// The guard initializes PETSc only when the embedding application has not
/// already done so.  It finalizes PETSc only when it owns that initialization.
///
/// Keep PetscContext alive longer than every PETSc object created by GSMakerLib.
class PetscContext final {
public:
    PetscContext();
    ~PetscContext() noexcept;

    PetscContext(
        const PetscContext&) = delete;

    PetscContext& operator=(
        const PetscContext&) = delete;

    PetscContext(
        PetscContext&&) = delete;

    PetscContext& operator=(
        PetscContext&&) = delete;

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


