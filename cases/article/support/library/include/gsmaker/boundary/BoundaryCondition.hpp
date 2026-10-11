#pragma once

#include <functional>
#include <variant>

#include "gsmaker/core/Types.hpp"


namespace gss {

/**
 * @brief Geometric point at which a boundary condition is evaluated.
 *
 * The point identifies the owner cell, the boundary-face position and the
 * outward unit normal.  The boundary-condition representation itself does
 * not depend on any mesh implementation.
 */
struct BoundaryPoint2D {
    Index owner_cell{};
    Real x{};
    Real y{};
    Real normal_x{};
    Real normal_y{};
};


/**
 * @brief Coefficients of the general boundary condition
 *
 *     alpha * phi + beta * d(phi)/dn = gamma.
 */
struct BoundaryCoefficients {
    Real alpha{};
    Real beta{};
    Real gamma{};
};


/**
 * @brief Spatially varying scalar boundary datum.
 *
 * std::function is used deliberately here because boundary data may be
 * constants, analytical functions, captured parameters or values supplied
 * by a study.  The public boundary-condition model remains value based and
 * uses no inheritance or virtual dispatch.
 */
using BoundaryDatum =
    std::function<Real(const BoundaryPoint2D&)>;


/**
 * @brief Dirichlet condition: phi = gamma.
 */
class Dirichlet {
public:
    explicit Dirichlet(
        BoundaryDatum gamma);

    explicit Dirichlet(
        Real gamma);

    [[nodiscard]] BoundaryCoefficients coefficients_at(
        const BoundaryPoint2D& point) const;

private:
    BoundaryDatum gamma_;
};


/**
 * @brief Neumann condition: d(phi)/dn = gamma.
 */
class Neumann {
public:
    explicit Neumann(
        BoundaryDatum gamma);

    explicit Neumann(
        Real gamma);

    [[nodiscard]] BoundaryCoefficients coefficients_at(
        const BoundaryPoint2D& point) const;

private:
    BoundaryDatum gamma_;
};


/**
 * @brief Robin condition:
 *
 *     alpha * phi + beta * d(phi)/dn = gamma.
 *
 * Both constant and spatially varying alpha/beta coefficients are supported.
 */
class Robin {
public:
    Robin(
        Real alpha,
        Real beta,
        BoundaryDatum gamma);

    Robin(
        Real alpha,
        Real beta,
        Real gamma);

    Robin(
        BoundaryDatum alpha,
        BoundaryDatum beta,
        BoundaryDatum gamma);

    [[nodiscard]] BoundaryCoefficients coefficients_at(
        const BoundaryPoint2D& point) const;

private:
    BoundaryDatum alpha_;
    BoundaryDatum beta_;
    BoundaryDatum gamma_;
};


/**
 * @brief Closed set of boundary-condition types used by GSMakerLib.
 *
 * The alternative carried by the variant is the condition type.  No base
 * class, virtual function or run-time type identification is required.
 */
using BoundaryCondition =
    std::variant<
        Dirichlet,
        Neumann,
        Robin>;


/**
 * @brief Evaluate any BoundaryCondition at a boundary point.
 */
[[nodiscard]] BoundaryCoefficients coefficients_at(
    const BoundaryCondition& condition,
    const BoundaryPoint2D& point);

} // namespace gss
