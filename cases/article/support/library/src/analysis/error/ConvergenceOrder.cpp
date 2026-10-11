#include "gsmaker/analysis/error/ConvergenceOrder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace gss {

namespace {

constexpr long double initial_order = 1.0L;
constexpr int maximum_iterations = 100;

constexpr long double function_tolerance = 1.0e-12L;
constexpr long double step_tolerance = 1.0e-12L;


[[nodiscard]] bool is_finite(
    Real value) noexcept {

    return std::isfinite(
        static_cast<double>(value));
}


[[nodiscard]] bool is_finite(
    long double value) noexcept {

    return std::isfinite(value);
}


void validate_result(
    Real value,
    const char* name) {

    if (!is_finite(value)) {
        throw std::invalid_argument(
            std::string{"ConvergenceOrder "} +
            name +
            " must be finite");
    }
}


[[nodiscard]] long double characteristic_size(
    Index nvol) {

    if (!(nvol > Index{0})) {
        throw std::invalid_argument(
            "ConvergenceOrder requires a positive number of finite volumes");
    }

    return 1.0L /
        std::sqrt(
            static_cast<long double>(
                nvol));
}


[[nodiscard]] long double equation_value(
    long double order,
    long double h1,
    long double h2,
    long double h3,
    long double result_ratio) {

    const long double h1p =
        std::pow(h1, order);

    const long double h2p =
        std::pow(h2, order);

    const long double h3p =
        std::pow(h3, order);

    const long double denominator =
        h3p - h2p;

    const long double denominator_scale =
        std::max(
            {
                std::abs(h2p),
                std::abs(h3p),
                1.0L
            });

    const long double denominator_limit =
        64.0L *
        std::numeric_limits<long double>::epsilon() *
        denominator_scale;

    if (std::abs(denominator) <= denominator_limit) {
        throw std::runtime_error(
            "ConvergenceOrder three-level equation is singular for the current order");
    }

    const long double value =
        (h1p - h2p) /
        denominator -
        result_ratio;

    if (!is_finite(value)) {
        throw std::runtime_error(
            "ConvergenceOrder produced a non-finite equation value");
    }

    return value;
}


[[nodiscard]] long double numerical_derivative(
    long double order,
    long double h1,
    long double h2,
    long double h3,
    long double result_ratio) {

    const long double epsilon =
        std::numeric_limits<long double>::epsilon();

    const long double delta =
        std::cbrt(epsilon) *
        std::max(
            1.0L,
            std::abs(order));

    const long double value_plus =
        equation_value(
            order + delta,
            h1,
            h2,
            h3,
            result_ratio);

    const long double value_minus =
        equation_value(
            order - delta,
            h1,
            h2,
            h3,
            result_ratio);

    const long double derivative =
        (value_plus - value_minus) /
        (2.0L * delta);

    if (!is_finite(derivative)) {
        throw std::runtime_error(
            "ConvergenceOrder produced a non-finite numerical derivative");
    }

    const long double derivative_limit =
        64.0L *
        std::sqrt(epsilon);

    if (std::abs(derivative) <= derivative_limit) {
        throw std::runtime_error(
            "ConvergenceOrder numerical derivative is too small for Newton-Raphson");
    }

    return derivative;
}

} // namespace


Real ConvergenceOrder::mesh_size(
    Index nvol) {

    return static_cast<Real>(
        characteristic_size(
            nvol));
}


Real ConvergenceOrder::estimate(
    Index nvol1,
    Real f1,
    Index nvol2,
    Real f2,
    Index nvol3,
    Real f3) {

    validate_result(
        f1,
        "f1");

    validate_result(
        f2,
        "f2");

    validate_result(
        f3,
        "f3");

    if (!(nvol1 > Index{0}) ||
        !(nvol2 > Index{0}) ||
        !(nvol3 > Index{0})) {

        throw std::invalid_argument(
            "ConvergenceOrder requires positive finite-volume counts");
    }

    if (!(nvol1 < nvol2 &&
          nvol2 < nvol3)) {

        throw std::invalid_argument(
            "ConvergenceOrder requires nvol1 < nvol2 < nvol3");
    }

    const long double h1 =
        characteristic_size(
            nvol1);

    const long double h2 =
        characteristic_size(
            nvol2);

    const long double h3 =
        characteristic_size(
            nvol3);

    const long double numerator =
        static_cast<long double>(f1) -
        static_cast<long double>(f2);

    const long double denominator =
        static_cast<long double>(f3) -
        static_cast<long double>(f2);

    const long double result_scale =
        std::max(
            {
                std::abs(
                    static_cast<long double>(f1)),
                std::abs(
                    static_cast<long double>(f2)),
                std::abs(
                    static_cast<long double>(f3)),
                1.0L
            });

    const long double result_limit =
        64.0L *
        std::numeric_limits<long double>::epsilon() *
        result_scale;

    if (std::abs(denominator) <= result_limit) {
        throw std::runtime_error(
            "ConvergenceOrder requires f3 - f2 to be non-zero");
    }

    const long double result_ratio =
        numerator /
        denominator;

    if (!is_finite(result_ratio)) {
        throw std::runtime_error(
            "ConvergenceOrder produced a non-finite result ratio");
    }

    long double order =
        initial_order;

    for (int iteration = 0;
         iteration < maximum_iterations;
         ++iteration) {

        const long double value =
            equation_value(
                order,
                h1,
                h2,
                h3,
                result_ratio);

        if (std::abs(value) <=
            function_tolerance) {

            return static_cast<Real>(
                order);
        }

        const long double derivative =
            numerical_derivative(
                order,
                h1,
                h2,
                h3,
                result_ratio);

        const long double correction =
            value /
            derivative;

        const long double next_order =
            order -
            correction;

        if (!is_finite(next_order)) {
            throw std::runtime_error(
                "ConvergenceOrder Newton-Raphson produced a non-finite order");
        }

        if (std::abs(correction) <=
            step_tolerance *
            std::max(
                1.0L,
                std::abs(next_order))) {

            const long double final_value =
                equation_value(
                    next_order,
                    h1,
                    h2,
                    h3,
                    result_ratio);

            if (std::abs(final_value) <=
                function_tolerance) {

                return static_cast<Real>(
                    next_order);
            }
        }

        order =
            next_order;
    }

    throw std::runtime_error(
        "ConvergenceOrder Newton-Raphson did not converge");
}

} // namespace gss
