#include "gsmaker/linear_system/PetscContext.hpp"

#include <sstream>
#include <stdexcept>
#include <string>

namespace gss {

void petsc_check(
    PetscErrorCode error,
    std::string_view action) {

    if (error == PETSC_SUCCESS) {
        return;
    }

    const char* text = nullptr;
    const char* specific = nullptr;

    (void)PetscErrorMessage(
        error,
        &text,
        nullptr); // Portable across PETSc versions; retain the generic error message.

    std::ostringstream message;

    message
        << "PETSc failure while attempting to "
        << action
        << " (error code "
        << error
        << ')';

    if (text != nullptr) {
        message
            << ": "
            << text;
    }

    if (specific != nullptr &&
        specific[0] != '\0') {

        message
            << " ["
            << specific
            << ']';
    }

    throw std::runtime_error(
        message.str());
}


PetscContext::PetscContext() {

    PetscBool finalized =
        PETSC_FALSE;

    petsc_check(
        PetscFinalized(
            &finalized),
        "query whether PETSc has been finalized");

    if (finalized == PETSC_TRUE) {
        throw std::runtime_error(
            "PETSc has already been finalized in this process");
    }

    PetscBool initialized_flag =
        PETSC_FALSE;

    petsc_check(
        PetscInitialized(
            &initialized_flag),
        "query whether PETSc is initialized");

    if (initialized_flag == PETSC_FALSE) {

        petsc_check(
            PetscInitializeNoArguments(),
            "initialize PETSc");

        owns_initialization_ =
            true;
    }
}


PetscContext::~PetscContext() noexcept {

    if (!owns_initialization_) {
        return;
    }

    PetscBool finalized =
        PETSC_FALSE;

    if (PetscFinalized(
            &finalized) !=
        PETSC_SUCCESS) {

        return;
    }

    if (finalized == PETSC_FALSE) {
        (void)PetscFinalize();
    }
}


bool PetscContext::initialized() {

    PetscBool initialized_flag =
        PETSC_FALSE;

    petsc_check(
        PetscInitialized(
            &initialized_flag),
        "query whether PETSc is initialized");

    return
        initialized_flag ==
        PETSC_TRUE;
}

} // namespace gss