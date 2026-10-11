find_package(MPI REQUIRED)
find_package(PkgConfig REQUIRED)

if(DEFINED ENV{PETSC_DIR})
  if(DEFINED ENV{PETSC_ARCH})
    set(ENV{PKG_CONFIG_PATH}
      "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
  endif()

  set(ENV{PKG_CONFIG_PATH}
    "$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()

file(GLOB GSS_LOCAL_PETSC_PKGCONFIG_DIRS
  "$ENV{HOME}/software/petsc/*/lib/pkgconfig")
foreach(petsc_pkgconfig_dir IN LISTS GSS_LOCAL_PETSC_PKGCONFIG_DIRS)
  set(ENV{PKG_CONFIG_PATH}
    "${petsc_pkgconfig_dir}:$ENV{PKG_CONFIG_PATH}")
endforeach()

pkg_check_modules(PETSC REQUIRED IMPORTED_TARGET PETSc)

if(NOT TARGET PETSc::PETSc)
  add_library(PETSc::PETSc INTERFACE IMPORTED)
  target_link_libraries(PETSc::PETSc INTERFACE PkgConfig::PETSC MPI::MPI_CXX)
endif()

find_package(ADOLC REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(OpenMP)
find_package(GTest REQUIRED)
find_package(VoronoiMeshMaker REQUIRED)
# Compatibility with the modular VoronoiMeshMaker package.
if(NOT TARGET VoronoiMeshMaker::VoronoiMeshMaker AND TARGET VoronoiMeshMaker::vmm)
  find_library(GSS_LEGACY_VMM_LIBRARY NAMES VoronoiMeshMaker)
  add_library(VoronoiMeshMaker::VoronoiMeshMaker INTERFACE IMPORTED)
  if(GSS_LEGACY_VMM_LIBRARY)
    # The legacy headers require the ErrorConfig implementation in this library.
    set_target_properties(VoronoiMeshMaker::VoronoiMeshMaker PROPERTIES
      INTERFACE_LINK_LIBRARIES "${GSS_LEGACY_VMM_LIBRARY};VoronoiMeshMaker::vmm")
  else()
    set_target_properties(VoronoiMeshMaker::VoronoiMeshMaker PROPERTIES
      INTERFACE_LINK_LIBRARIES VoronoiMeshMaker::vmm)
  endif()
endif()
find_package(VoroGradToolkit REQUIRED)
find_package(VoronoiGradientKit QUIET)

message(STATUS "MPI CXX compiler   : ${MPI_CXX_COMPILER}")
message(STATUS "PETSc version      : ${PETSC_VERSION}")
message(STATUS "VoronoiMeshMaker   : ${VoronoiMeshMaker_FOUND}")
message(STATUS "VoroGradToolkit    : ${VoroGradToolkit_FOUND}")
message(STATUS "VoronoiGradientKit : ${VoronoiGradientKit_FOUND}")
