
#pragma once

#include "gsmaker/core/Types.hpp"

namespace gss {

/**
 * @brief Observed convergence order obtained from three mesh levels.
 *
 * For three calculations (N1,f1), (N2,f2), (N3,f3), the characteristic
 * mesh sizes are defined by
 *
 *     h_i = 1 / sqrt(N_i),
 *
 * where N_i is the number of finite volumes.
 *
 * Assuming
 *
 *     f(h) = f_exact + C h^p,
 *
 * the observed order p satisfies
 *
 *     (h1^p - h2^p) / (h3^p - h2^p)
 *         =
 *     (f1 - f2) / (f3 - f2).
 *
 * The nonlinear scalar equation is solved with Newton-Raphson, starting from
 * p0 = 1.  The derivative with respect to p is evaluated numerically by a
 * centred finite difference.
 */
class ConvergenceOrder {
public:
    /**
     * @brief Estimate the observed order from three successively refined meshes.
     *
     * The meshes must satisfy
     *
     *     nvol1 < nvol2 < nvol3,
     *
     * so that
     *
     *     h1 > h2 > h3.
     *
     * @param nvol1 Number of finite volumes in the coarsest mesh.
     * @param f1    Scalar result obtained on the coarsest mesh.
     * @param nvol2 Number of finite volumes in the intermediate mesh.
     * @param f2    Scalar result obtained on the intermediate mesh.
     * @param nvol3 Number of finite volumes in the finest mesh.
     * @param f3    Scalar result obtained on the finest mesh.
     *
     * @return Observed convergence order p.
     *
     * @throws std::invalid_argument for invalid mesh sizes or non-finite data.
     * @throws std::runtime_error if the three-level equation is degenerate or
     * Newton-Raphson does not converge.
     */
    [[nodiscard]] static Real estimate(
        Index nvol1,
        Real f1,
        Index nvol2,
        Real f2,
        Index nvol3,
        Real f3);

    /**
     * @brief Characteristic mesh size used by the convergence analysis.
     *
     *     h = 1 / sqrt(Nvol)
     */
    [[nodiscard]] static Real mesh_size(
        Index nvol);
};

} // namespace gss
