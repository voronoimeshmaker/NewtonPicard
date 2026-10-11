\
    #pragma once

    #include <concepts>

    #include "gsmaker/core/Types.hpp"

    namespace gss {

    /**
     * @brief C++20 contract for a two-dimensional analytical reference solution.
     *
     * An analytical solution used by the verification layer must provide the
     * scalar field psi(x,y).  Problem-specific auxiliary functions, such as the
     * Grad-Shafranov free function F(Psi), deliberately remain outside this
     * minimal contract.
     *
     * The concept is intentionally non-intrusive: analytical reference
     * solutions do not participate in the solver hierarchy and require no
     * run-time polymorphism.
     */
    template<typename T>
    concept AnalyticalSolution =
        requires(const T& solution, Real x, Real y) {
            {
                solution.psi(x, y)
            } -> std::convertible_to<Real>;
        };

    } // namespace gss

    