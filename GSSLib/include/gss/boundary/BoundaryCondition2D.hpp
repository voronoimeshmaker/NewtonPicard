// SPDX-License-Identifier: GPL-3.0-or-later
/** @file BoundaryCondition2D.hpp
 * @brief Value-based two-dimensional boundary conditions without inheritance.
 * @ingroup boundary
 */
#pragma once

#include <concepts>
#include <functional>
#include <utility>

#include "gss/core/Types.hpp"

namespace gss {

class BoundaryApplicationPoint2D {
public:
    BoundaryApplicationPoint2D(Index volume_id, Real x, Real y,
                               Real outward_normal_x,
                               Real outward_normal_y) noexcept
        : volume_id_(volume_id), x_(x), y_(y),
          outward_normal_x_(outward_normal_x),
          outward_normal_y_(outward_normal_y) {}

    [[nodiscard]] Index volume_id() const noexcept { return volume_id_; }
    [[nodiscard]] Real x() const noexcept { return x_; }
    [[nodiscard]] Real y() const noexcept { return y_; }
    [[nodiscard]] Real outward_normal_x() const noexcept { return outward_normal_x_; }
    [[nodiscard]] Real outward_normal_y() const noexcept { return outward_normal_y_; }

private:
    Index volume_id_{};
    Real x_{};
    Real y_{};
    Real outward_normal_x_{};
    Real outward_normal_y_{};
};

struct BoundaryCoefficients2D {
    Real alpha{};
    Real beta{};
    Real gamma{};
};

template<class Condition>
concept BoundaryCondition2D = requires(
    const Condition& condition,
    const BoundaryApplicationPoint2D& point) {
    { condition.coefficients_at(point) }
        -> std::same_as<BoundaryCoefficients2D>;
};

using BoundaryDatum2D = std::function<Real(const BoundaryApplicationPoint2D&)>;

class DirichletBoundary2D {
public:
    explicit DirichletBoundary2D(BoundaryDatum2D gamma)
        : gamma_(std::move(gamma)) {}

    [[nodiscard]] BoundaryCoefficients2D coefficients_at(
        const BoundaryApplicationPoint2D& point) const {
        return {Real{1}, Real{0}, gamma_(point)};
    }

private:
    BoundaryDatum2D gamma_;
};

class NeumannBoundary2D {
public:
    explicit NeumannBoundary2D(BoundaryDatum2D gamma)
        : gamma_(std::move(gamma)) {}

    [[nodiscard]] BoundaryCoefficients2D coefficients_at(
        const BoundaryApplicationPoint2D& point) const {
        return {Real{0}, Real{1}, gamma_(point)};
    }

private:
    BoundaryDatum2D gamma_;
};

class RobinBoundary2D {
public:
    RobinBoundary2D(Real alpha, Real beta, BoundaryDatum2D gamma)
        : alpha_(alpha), beta_(beta), gamma_(std::move(gamma)) {}

    [[nodiscard]] BoundaryCoefficients2D coefficients_at(
        const BoundaryApplicationPoint2D& point) const {
        return {alpha_, beta_, gamma_(point)};
    }

private:
    Real alpha_{};
    Real beta_{};
    BoundaryDatum2D gamma_;
};

static_assert(BoundaryCondition2D<DirichletBoundary2D>);
static_assert(BoundaryCondition2D<NeumannBoundary2D>);
static_assert(BoundaryCondition2D<RobinBoundary2D>);

} // namespace gss
