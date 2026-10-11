#include "gsmaker/analysis/error/ErrorNorms.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace gss {

namespace {

[[nodiscard]] bool is_finite(Real value) noexcept {

    return std::isfinite(
        static_cast<double>(value));
}


void validate_error(
    std::span<const Real> error) {

    if (error.empty()) {
        throw std::invalid_argument(
            "ErrorNorms requires a non-empty error field");
    }

    for (const Real value : error) {

        if (!is_finite(value)) {
            throw std::invalid_argument(
                "ErrorNorms error field contains a non-finite value");
        }
    }
}


void validate_inputs(
    std::span<const Real> error,
    std::span<const Real> cell_volume) {

    validate_error(error);

    if (cell_volume.size() != error.size()) {
        throw std::invalid_argument(
            "ErrorNorms error and cell-volume arrays must have the same size");
    }

    for (const Real volume : cell_volume) {

        if (!is_finite(volume)) {
            throw std::invalid_argument(
                "ErrorNorms cell-volume array contains a non-finite value");
        }

        if (!(volume > Real{0})) {
            throw std::invalid_argument(
                "ErrorNorms requires strictly positive cell volumes");
        }
    }
}

} // namespace


ErrorNorms::Values ErrorNorms::compute(
    std::span<const Real> error,
    std::span<const Real> cell_volume) {

    validate_inputs(
        error,
        cell_volume);

    long double weighted_l1 = 0.0L;
    long double weighted_l2_squared = 0.0L;
    long double total_volume = 0.0L;
    long double maximum_absolute_error = 0.0L;

    for (std::size_t i = 0; i < error.size(); ++i) {

        const long double e =
            static_cast<long double>(
                error[i]);

        const long double volume =
            static_cast<long double>(
                cell_volume[i]);

        const long double absolute_error =
            std::abs(e);

        total_volume += volume;

        weighted_l1 +=
            volume * absolute_error;

        weighted_l2_squared +=
            volume * e * e;

        maximum_absolute_error =
            std::max(
                maximum_absolute_error,
                absolute_error);
    }

    const long double weighted_l2 =
        std::sqrt(
            weighted_l2_squared);

    const long double weighted_rms =
        std::sqrt(
            weighted_l2_squared /
            total_volume);

    return Values{
        static_cast<Real>(weighted_l1),
        static_cast<Real>(weighted_l2),
        static_cast<Real>(maximum_absolute_error),
        static_cast<Real>(weighted_rms),
        static_cast<Real>(total_volume)
    };
}


Real ErrorNorms::l1(
    std::span<const Real> error,
    std::span<const Real> cell_volume) {

    return compute(
        error,
        cell_volume).l1;
}


Real ErrorNorms::l2(
    std::span<const Real> error,
    std::span<const Real> cell_volume) {

    return compute(
        error,
        cell_volume).l2;
}


Real ErrorNorms::linf(
    std::span<const Real> error) {

    validate_error(error);

    long double maximum_absolute_error = 0.0L;

    for (const Real value : error) {

        maximum_absolute_error =
            std::max(
                maximum_absolute_error,
                std::abs(
                    static_cast<long double>(
                        value)));
    }

    return static_cast<Real>(
        maximum_absolute_error);
}


Real ErrorNorms::rms(
    std::span<const Real> error,
    std::span<const Real> cell_volume) {

    return compute(
        error,
        cell_volume).rms;
}

} // namespace gss
