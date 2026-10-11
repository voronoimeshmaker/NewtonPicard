\
#include "gsmaker/analysis/error/ReferenceError.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace gss {

namespace {

[[nodiscard]] bool is_finite(Real value) noexcept {

    return std::isfinite(
        static_cast<double>(value));
}


void validate_field(
    std::span<const Real> field,
    const char* field_name) {

    for (const Real value : field) {

        if (!is_finite(value)) {
            throw std::invalid_argument(
                std::string{"ReferenceError "} +
                field_name +
                " field contains a non-finite value");
        }
    }
}

} // namespace


std::vector<Real> ReferenceError::compute(
    std::span<const Real> numerical,
    std::span<const Real> reference) {

    if (numerical.empty()) {
        throw std::invalid_argument(
            "ReferenceError requires non-empty fields");
    }

    if (numerical.size() != reference.size()) {
        throw std::invalid_argument(
            "ReferenceError numerical and reference fields must have the same size");
    }

    validate_field(
        numerical,
        "numerical");

    validate_field(
        reference,
        "reference");

    std::vector<Real> error(
        numerical.size());

    for (std::size_t i = 0; i < numerical.size(); ++i) {

        error[i] =
            numerical[i] -
            reference[i];
    }

    return error;
}

} // namespace gss
