#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <VoroGradToolkit/Fields/ScalarField2D.hpp>
#include <VoroGradToolkit/Gradients/Common/GradientRunResult2D.hpp>
#include <VoroGradToolkit/Mesh/GradientWorkspace2D.hpp>

#include "gsmaker/core/Types.hpp"

namespace gss {

    
/// @brief Solver-facing gradient field in finite-volume (volume-id) order.
struct GradientField2D {
    std::vector<Real> x;
    std::vector<Real> y;

    [[nodiscard]]
    std::size_t size() const noexcept {
        return x.size();
    }
};

/// @brief Diagnostics returned by the external reconstruction method.
struct GradientDiagnostics2D {
    std::size_t processed_volumes{};
    std::size_t failed_volumes{};
    std::size_t singular_volumes{};
};

struct GradientResult2D {
    GradientField2D field;
    GradientDiagnostics2D diagnostics;
};

namespace detail {

template<typename Method>
[[nodiscard]]
auto call_reconstruct_gradients(
    const vgt::GradientWorkspace2D& workspace,
    const vgt::ScalarField2D& scalar,
    const Method& method,
    vgt::GradientField2D& output)
    -> decltype(
        reconstruct_gradients(
            workspace,
            scalar,
            method,
            output))
{
    // Deliberately unqualified.  Argument-dependent lookup selects the
    // reconstruction supplied by VoroGradToolkit for Method.
    return reconstruct_gradients(
        workspace,
        scalar,
        method,
        output);
}

template<typename Method>
concept ToolkitGradientMethod2D =
    requires(
        const vgt::GradientWorkspace2D& workspace,
        const vgt::ScalarField2D& scalar,
        const Method& method,
        vgt::GradientField2D& output) {

        {
            call_reconstruct_gradients(
                workspace,
                scalar,
                method,
                output)
        } -> std::same_as<vgt::gradients::GradientRunResult2D>;
    };

[[nodiscard]]
inline bool finite(
    Real value) noexcept {

    return std::isfinite(
        static_cast<double>(value));
}

} // namespace detail


/// @brief Generic bridge from GSMaker finite-volume fields to VoroGradToolkit.
///
/// Gradient contains no reconstruction algorithm and no method enumeration.
/// The concrete method is supplied by VoroGradToolkit and selected at the
/// call site through C++20 static polymorphism.
///
/// Any toolkit object for which ADL can resolve
///
///     reconstruct_gradients(workspace, scalar, method, output)
///
/// is accepted.  This includes, for example, WLSQR2D, WLSQRCache2D and
/// ModifiedGreenGauss2D without any switch/enum/inheritance inside GSMakerLib.
class Gradient final {
public:
    template<detail::ToolkitGradientMethod2D Method>
    [[nodiscard]]
    static GradientResult2D compute(
        const vgt::GradientWorkspace2D& workspace,
        std::span<const Real> field,
        const Method& method) {

        const std::size_t n =
            workspace.volume_count();

        if (n == 0U ||
            field.size() != n) {

            throw std::invalid_argument(
                "Gradient field size must equal the non-zero workspace volume count");
        }

        for (const Real value : field) {
            if (!detail::finite(value)) {
                throw std::invalid_argument(
                    "Gradient input field contains a non-finite value");
            }
        }

        // VoroGradToolkit scalar fields are indexed by site.  GSMakerLib
        // scientific fields are indexed by finite-volume id, so convert
        // explicitly using the mapping already owned by GradientWorkspace2D.
        vgt::ScalarField2D scalar(
            workspace.site_count());

        for (const auto volume_id :
             workspace.all_volume_ids()) {

            const auto& volume =
                workspace.volume(
                    volume_id);

            scalar[volume.site_index] =
                field[volume_id];
        }

        vgt::GradientField2D external(
            workspace.volume_count());

        const auto run =
            detail::call_reconstruct_gradients(
                workspace,
                scalar,
                method,
                external);

        GradientResult2D result{
            GradientField2D{
                std::vector<Real>(n, Real{0}),
                std::vector<Real>(n, Real{0})},
            GradientDiagnostics2D{
                run.diagnostics.processed_volumes,
                run.diagnostics.failed_volumes,
                run.diagnostics.singular_volumes}
        };

        for (const auto volume_id :
             workspace.all_volume_ids()) {

            const auto& value =
                external[volume_id];

            if (!detail::finite(value.x) ||
                !detail::finite(value.y)) {

                throw std::runtime_error(
                    "Gradient reconstruction produced a non-finite value");
            }

            result.field.x[volume_id] =
                static_cast<Real>(
                    value.x);

            result.field.y[volume_id] =
                static_cast<Real>(
                    value.y);
        }

        return result;
    }
};

static_assert(
    !std::is_polymorphic_v<Gradient>);

} // namespace gss
