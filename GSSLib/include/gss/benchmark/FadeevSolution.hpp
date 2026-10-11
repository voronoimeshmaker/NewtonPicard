#pragma once

#include <cmath>
#include <numbers>
#include <utility>
#include <vector>

#include "gss/core/Types.hpp"

namespace gss {

/// @brief Fadeev analytical magnetic-island solution.
struct FadeevSolution {
    Real epsilon{0.1};
    Real psi0{1.0};
    Real k{1.0};

    /// @brief Exact scalar potential.
    [[nodiscard]] Real psi_exact(Real x, Real y) const noexcept {
        return psi0 * std::log(std::cosh(y) - epsilon * std::cos(k * x)) - std::log(1.0 - epsilon);
    }

    /// @brief x derivative of exact scalar potential.
    [[nodiscard]] Real dpsi_dx(Real x, Real y) const noexcept {
        const Real q = std::cosh(y) - epsilon * std::cos(k * x);
        return psi0 * epsilon * k * std::sin(k * x) / q;
    }

    /// @brief y derivative of exact scalar potential.
    [[nodiscard]] Real dpsi_dy(Real x, Real y) const noexcept {
        const Real q = std::cosh(y) - epsilon * std::cos(k * x);
        return psi0 * std::sinh(y) / q;
    }

    /// @brief xx derivative of exact scalar potential.
    [[nodiscard]] Real d2psi_dx2(Real x, Real y) const noexcept {
        const Real c = std::cos(k * x);
        const Real s = std::sin(k * x);
        const Real q = std::cosh(y) - epsilon * c;
        return psi0 * epsilon * k * k * (c * q - epsilon * s * s) / (q * q);
    }

    /// @brief yy derivative of exact scalar potential.
    [[nodiscard]] Real d2psi_dy2(Real x, Real y) const noexcept {
        const Real sh = std::sinh(y);
        const Real ch = std::cosh(y);
        const Real q = ch - epsilon * std::cos(k * x);
        return psi0 * (ch * q - sh * sh) / (q * q);
    }

    /// @brief Exact Laplacian.
    [[nodiscard]] Real laplacian(Real x, Real y) const noexcept {
        return d2psi_dx2(x, y) + d2psi_dy2(x, y);
    }

    /// @brief Exact nonlinear source term.
    [[nodiscard]] Real source_exact(Real x, Real y) const noexcept {
        return laplacian(x, y);
    }

    /// @brief Analytical X-point locations in the default first periodic cell.
    [[nodiscard]] std::vector<std::pair<Real, Real>> xpoint_locations() const {
        return {{0.0, 0.0}, {2.0 * std::numbers::pi_v<Real> / k, 0.0}};
    }

    /// @brief Analytical O-point location in the default first periodic cell.
    [[nodiscard]] std::vector<std::pair<Real, Real>> opoint_locations() const {
        return {{std::numbers::pi_v<Real> / k, 0.0}};
    }
};

} // namespace gss
