find_path(PETSC_INCLUDE_DIR petsc.h
  HINTS "$ENV{PETSC_DIR}/include" "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include")

find_library(PETSC_LIBRARY petsc
  HINTS "$ENV{PETSC_DIR}/lib" "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc REQUIRED_VARS PETSC_INCLUDE_DIR PETSC_LIBRARY)

if(PETSc_FOUND AND NOT TARGET PETSc::PETSc)
  add_library(PETSc::PETSc UNKNOWN IMPORTED)
  set_target_properties(PETSc::PETSc PROPERTIES
    IMPORTED_LOCATION "${PETSC_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIR}")
endif()
