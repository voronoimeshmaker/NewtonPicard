#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

#include "gsmaker/core/Types.hpp"
#include "gsmaker/gradient/Gradient.hpp"

namespace gss {

struct GradientRatioSourceResult {
    std::vector<Real> values;
    GradientField2D grad_psi;
    GradientField2D grad_F;

    
    std::size_t regularised_cells{};

    /// One entry per cell, non-zero where the threshold altered the value.
    /// Those entries are not recovered derivatives and must be either
    /// discarded or repaired before use; see StationaryPointRecovery.
    std::vector<std::uint8_t> regularised;

    /// Largest |grad Psi|^2 over the mesh; the scale the threshold is
    /// measured against.
    Real gradient_scale{};

    /// Threshold actually applied, equal to regularisation_delta^2 times
    /// gradient_scale.
    Real regularisation_threshold{};

    /// Smallest |grad Psi|^2 relative to gradient_scale.  A value approaching
    /// zero means the mesh resolves a stationary point of Psi.
    Real smallest_relative_gradient{};

    GradientDiagnostics2D psi_gradient_diagnostics{};
    GradientDiagnostics2D free_function_gradient_diagnostics{};
};

/// @brief Recover dF/dPsi numerically from externally reconstructed gradients.
///
/// The gradient reconstruction method is supplied by VoroGradToolkit:
///
///     dF/dPsi = (grad F . grad Psi) / (|grad Psi|^2 + tau^2)
///
/// GradientRatioSource deliberately knows nothing about WLSQR, MGG or any
/// other concrete reconstruction algorithm.
///
/// The regularisation parameter is RELATIVE, not absolute: the threshold
/// applied is
///
///     tau^2 = regularisation_delta^2 * max_i |grad Psi_i|^2
///
/// The recovered derivative is a ratio of two gradients and therefore carries
/// no units of Psi.  An absolute threshold destroys that invariance: it leaves
/// the ratio untouched while |grad Psi|^2 falls below it, so the same physical
/// problem stated in different units returns a different, silently plausible
/// answer.
///
/// CONDITIONING.  The ratio is well defined only where grad Psi does not
/// vanish.  With gradients reconstructed to second order, the recovered
/// derivative carries a relative error of order h^2 / |grad Psi|: it is second
/// order where the flux gradient is O(1), degrades to first order in the
/// neighbourhood of a stationary point, where |grad Psi| = O(h), and is
/// undefined at the stationary point itself, where the ratio is 0/0.  No
/// choice of threshold removes that indeterminacy; it only decides whether the
/// caller receives a zero, a noise-dominated value, or a flag.  This class
/// raises the flag.
class GradientRatioSource final {
public:
    template<detail::ToolkitGradientMethod2D Method>
    [[nodiscard]]
    static GradientRatioSourceResult compute(
        const vgt::GradientWorkspace2D& workspace,
        std::span<const Real> psi,
        std::span<const Real> free_function,
        const Method& method,
        Real regularisation_delta) {

        const std::size_t n =
            workspace.volume_count();

        if (n == 0U ||
            psi.size() != n ||
            free_function.size() != n) {

            throw std::invalid_argument(
                "GradientRatioSource arrays must match the non-zero workspace volume count");
        }

        if (!std::isfinite(
                static_cast<double>(
                    regularisation_delta)) ||
            regularisation_delta < Real{0}) {

            throw std::invalid_argument(
                "GradientRatioSource regularisation delta must be finite and non-negative");
        }

        const auto psi_gradient =
            Gradient::compute(
                workspace,
                psi,
                method);

        const auto free_function_gradient =
            Gradient::compute(
                workspace,
                free_function,
                method);

        if (psi_gradient.diagnostics.failed_volumes != 0U ||
            free_function_gradient.diagnostics.failed_volumes != 0U) {

            throw std::runtime_error(
                "GradientRatioSource received a failed gradient reconstruction from VoroGradToolkit");
        }

        GradientRatioSourceResult result;
        result.grad_psi =
            psi_gradient.field;

        result.grad_F =
            free_function_gradient.field;

        result.psi_gradient_diagnostics =
            psi_gradient.diagnostics;

        result.free_function_gradient_diagnostics =
            free_function_gradient.diagnostics;

        result.values.resize(
            n,
            Real{0});

        result.regularised.assign(
            n,
            std::uint8_t{0});

        // First pass: the scale of the flux gradient over this mesh.  The
        // threshold is meaningless before it is known, which is why it cannot
        // be a constant supplied by configuration.
        Real largest = Real{0};
        Real smallest =
            std::numeric_limits<Real>::max();

        for (std::size_t i = 0;
             i < n;
             ++i) {

            const Real gx =
                result.grad_psi.x[i];

            const Real gy =
                result.grad_psi.y[i];

            const Real g2 =
                gx * gx +
                gy * gy;

            largest = std::max(largest, g2);
            smallest = std::min(smallest, g2);
        }

        result.gradient_scale = largest;

        result.regularisation_threshold =
            regularisation_delta *
            regularisation_delta *
            largest;

        result.smallest_relative_gradient =
            largest > Real{0}
                ? smallest / largest
                : Real{0};

        const Real threshold =
            result.regularisation_threshold;

        for (std::size_t i = 0;
             i < n;
             ++i) {

            const Real gx =
                result.grad_psi.x[i];

            const Real gy =
                result.grad_psi.y[i];

            const Real g2 =
                gx * gx +
                gy * gy;

            const Real numerator =
                result.grad_F.x[i] * gx +
                result.grad_F.y[i] * gy;

            const Real denominator =
                g2 + threshold;

            if (g2 <= threshold) {
                ++result.regularised_cells;
                result.regularised[i] = std::uint8_t{1};
            }

            result.values[i] =
                denominator == Real{0}
                    ? Real{0}
                    : numerator / denominator;

            if (!std::isfinite(
                    static_cast<double>(
                        result.values[i]))) {

                throw std::runtime_error(
                    "GradientRatioSource produced a non-finite source derivative");
            }
        }

        return result;
    }
};

} // namespace gss