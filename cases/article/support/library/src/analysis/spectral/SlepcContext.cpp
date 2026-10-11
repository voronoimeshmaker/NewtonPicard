#ifdef GSMAKER_USE_SLEPC

#include "gsmaker/analysis/spectral/SlepcContext.hpp"

#include <stdexcept>

#include "gsmaker/linear_system/PetscContext.hpp"

namespace gss {

    
SlepcContext::SlepcContext() {

    PetscBool finalized = PETSC_FALSE;

    petsc_check(
        SlepcFinalized(&finalized),
        "query whether SLEPc has been finalized");

    if (finalized == PETSC_TRUE) {
        throw std::runtime_error(
            "SLEPc has already been finalized in this process");
    }

    PetscBool initialized_flag = PETSC_FALSE;

    petsc_check(
        SlepcInitialized(&initialized_flag),
        "query whether SLEPc is initialized");

    if (initialized_flag == PETSC_FALSE) {
        petsc_check(
            SlepcInitializeNoArguments(),
            "initialize SLEPc");

        owns_initialization_ = true;
    }
}


SlepcContext::~SlepcContext() noexcept {

    if (!owns_initialization_) {
        return;
    }

    PetscBool finalized = PETSC_FALSE;

    if (SlepcFinalized(&finalized) != PETSC_SUCCESS) {
        return;
    }

    if (finalized == PETSC_FALSE) {
        (void)SlepcFinalize();
    }
}


bool SlepcContext::initialized() {

    PetscBool initialized_flag = PETSC_FALSE;

    petsc_check(
        SlepcInitialized(&initialized_flag),
        "query whether SLEPc is initialized");

    return initialized_flag == PETSC_TRUE;
}

} // namespace gss

#endif // GSMAKER_USE_SLEPC
