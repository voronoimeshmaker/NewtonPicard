#pragma once

#include <string_view>

#include "gss/solvers/Tags.hpp"

namespace gss {

template<typename Tag>
struct LinearSolverTraits;

template<>
struct LinearSolverTraits<CholeskySparseTag> {
    static constexpr std::string_view name = "PETSc sparse Cholesky";
    static constexpr bool direct = true;
};

template<>
struct LinearSolverTraits<AMGTag> {
    static constexpr std::string_view name = "PETSc GAMG";
    static constexpr bool direct = false;
};

// EXTENSION POINT: add new Tag specialization here
// See docs/extending.rst for instructions

} // namespace gss
