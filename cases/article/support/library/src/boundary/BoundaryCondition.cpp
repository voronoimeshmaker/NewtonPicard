#include "gsmaker/boundary/BoundaryCondition.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace gss {

namespace {

[[nodiscard]] bool is_finite(
    Real value) noexcept {

    return std::isfinite(
        static_cast<double>(value));
}


void require_datum(
    const BoundaryDatum& datum,
    const char* name) {

    if (!datum) {
        throw std::invalid_argument(
            std::string{"BoundaryCondition requires a valid "} +
            name +
            " datum");
    }
}



void validate_coefficients(
    const BoundaryCoefficients& coefficients) {

    if (!is_finite(coefficients.alpha) ||
        !is_finite(coefficients.beta) ||
        !is_finite(coefficients.gamma)) {

        throw std::runtime_error(
            "BoundaryCondition evaluation produced a non-finite coefficient");
    }

    if (coefficients.alpha == Real{0} &&
        coefficients.beta == Real{0}) {

        throw std::runtime_error(
            "BoundaryCondition requires alpha and beta not to be zero simultaneously");
    }
}


[[nodiscard]] BoundaryDatum constant_datum(
    Real value) {

    if (!is_finite(value)) {
        throw std::invalid_argument(
            "BoundaryCondition constant datum must be finite");
    }

    return [value](
        const BoundaryPoint2D&) noexcept {

        return value;
    };
}

} // namespace


Dirichlet::Dirichlet(
    BoundaryDatum gamma)
    : gamma_(std::move(gamma)) {

    require_datum(
        gamma_,
        "Dirichlet gamma");
}


Dirichlet::Dirichlet(
    Real gamma)
    : Dirichlet(
          constant_datum(gamma)) {}


BoundaryCoefficients Dirichlet::coefficients_at(
    const BoundaryPoint2D& point) const {

    const BoundaryCoefficients coefficients{
        Real{1},
        Real{0},
        gamma_(point)
    };

    validate_coefficients(
        coefficients);

    return coefficients;
}


Neumann::Neumann(
    BoundaryDatum gamma)
    : gamma_(std::move(gamma)) {

    require_datum(
        gamma_,
        "Neumann gamma");
}


Neumann::Neumann(
    Real gamma)
    : Neumann(
          constant_datum(gamma)) {}


BoundaryCoefficients Neumann::coefficients_at(
    const BoundaryPoint2D& point) const {

    const BoundaryCoefficients coefficients{
        Real{0},
        Real{1},
        gamma_(point)
    };

    validate_coefficients(
        coefficients);

    return coefficients;
}


Robin::Robin(
    Real alpha,
    Real beta,
    BoundaryDatum gamma)
    : Robin(
          constant_datum(alpha),
          constant_datum(beta),
          std::move(gamma)) {}


Robin::Robin(
    Real alpha,
    Real beta,
    Real gamma)
    : Robin(
          constant_datum(alpha),
          constant_datum(beta),
          constant_datum(gamma)) {}


Robin::Robin(
    BoundaryDatum alpha,
    BoundaryDatum beta,
    BoundaryDatum gamma)
    : alpha_(std::move(alpha)),
      beta_(std::move(beta)),
      gamma_(std::move(gamma)) {

    require_datum(
        alpha_,
        "Robin alpha");

    require_datum(
        beta_,
        "Robin beta");

    require_datum(
        gamma_,
        "Robin gamma");
}


BoundaryCoefficients Robin::coefficients_at(
    const BoundaryPoint2D& point) const {

    const BoundaryCoefficients coefficients{
        alpha_(point),
        beta_(point),
        gamma_(point)
    };

    validate_coefficients(
        coefficients);

    return coefficients;
}


BoundaryCoefficients coefficients_at(
    const BoundaryCondition& condition,
    const BoundaryPoint2D& point) {

    return std::visit(
        [&point](const auto& concrete_condition) {
            return concrete_condition.coefficients_at(
                point);
        },
        condition);
}

} // namespace gss
