#pragma once

#include <string_view>

#include "gss/solvers/Tags.hpp"

namespace gss {

template<typename Tag>
struct NonlinearSolverTraits;

template<>
struct NonlinearSolverTraits<PicardTag> {
    static constexpr bool needs_jacobian = false;
    static constexpr bool matrix_is_fixed = true;
    static constexpr std::string_view name = "Picard";
};

template<>
struct NonlinearSolverTraits<NewtonTag> {
    static constexpr bool needs_jacobian = true;
    static constexpr bool matrix_is_fixed = false;
    static constexpr std::string_view name = "Newton";
};

// EXTENSION POINT: add new Tag specialization here
// See docs/extending.rst for instructions

template<typename Backend>
concept LinearSolverConcept = requires(Backend backend, Mat A, Vec b, Vec x) {
    { backend.factorize(A) };
    { backend.solve(b, x) };
};

} // namespace gss
