cmake_minimum_required(VERSION 3.20)

if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
  set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
endif()

project(GradShafranovSolver VERSION 0.1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

include(Dependencies)
include(CompilerFlags)

file(GLOB_RECURSE GSS_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/GSSLib/include/*.hpp")
file(GLOB_RECURSE GSS_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/GSSLib/src/*.cpp")
list(FILTER GSS_SOURCES EXCLUDE REGEX ".*/main\\.cpp$")

add_library(gss ${GSS_HEADERS} ${GSS_SOURCES})
add_library(gss::gss ALIAS gss)
target_include_directories(gss
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GSSLib/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(gss
  PUBLIC PETSc::PETSc yaml-cpp::yaml-cpp ADOLC::ADOLC
  PUBLIC
    VoronoiMeshMaker::VoronoiMeshMaker
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    $<$<BOOL:${VoronoiGradientKit_FOUND}>:VoronoiGradientKit::VoronoiGradientKit>
  PRIVATE
    $<$<BOOL:${CGAL_FOUND}>:CGAL::CGAL>)
gss_apply_common_flags(gss)

if(CGAL_FOUND)
  target_compile_definitions(gss PUBLIC USE_CGAL=1)
endif()

if(VoronoiGradientKit_FOUND)
  target_compile_definitions(gss PUBLIC GSS_USE_VORONOI_GRADIENT_KIT=1)
endif()

add_executable(gss_solver GSSLib/src/main.cpp)
target_link_libraries(gss_solver PRIVATE gss::gss)
gss_apply_common_flags(gss_solver)

file(GLOB EXAMPLE_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cpp")
foreach(example_source IN LISTS EXAMPLE_SOURCES)
  get_filename_component(example_name "${example_source}" NAME_WE)
  add_executable("${example_name}" "${example_source}")
  target_link_libraries("${example_name}" PRIVATE gss::gss)
  gss_apply_common_flags("${example_name}")
endforeach()

file(GLOB_RECURSE PAPER_PROGRAM_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/programs/paper/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/programs/paper_analytic_source_test/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/programs/batch3_boundary_condition_20260903/src/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/programs/batch4_article_release_20260904/src/*.cpp")
foreach(program_source IN LISTS PAPER_PROGRAM_SOURCES)
  get_filename_component(program_name "${program_source}" NAME_WE)
  add_executable("${program_name}" "${program_source}")
  target_link_libraries("${program_name}" PRIVATE gss::gss)
  gss_apply_common_flags("${program_name}")
endforeach()

include(CTest)
if(BUILD_TESTING AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
  add_subdirectory(tests)
endif()

install(TARGETS gss gss_solver
  EXPORT GradShafranovSolverTargets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/GSSLib/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING
    PATTERN "*.hpp"
    PATTERN "*.h"
    PATTERN "*.tpp")

install(EXPORT GradShafranovSolverTargets
  NAMESPACE GradShafranovSolver::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GradShafranovSolver")

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/GradShafranovSolverConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GradShafranovSolverConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GradShafranovSolver")

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/GradShafranovSolverConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)

install(FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindADOLC.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/GradShafranovSolverConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/GradShafranovSolverConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GradShafranovSolver")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/verification/CMakeLists.txt")
  add_subdirectory(verification)
endif()

add_subdirectory(cases/fadeev)
