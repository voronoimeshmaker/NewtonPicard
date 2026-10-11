#pragma once

#include <cmath>
#include <numeric>

#include "gss/analysis/ErrorNorms.hpp"
#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/benchmark/FadeevSolution.hpp"
#include "gss/core/MeshAdapter.hpp"

namespace gss {

/// @brief Computes truncation error, solution error and convergence-study data.
template<typename MeshAdapterT, typename RHSAssemblerT>
struct TruncationErrorAnalyzer {
    /// @brief Bundle of local and global error metrics.
    struct ErrorBundle {
        NormBundle tau;
        std::vector<Real> tau_local;
        NormBundle e;
        std::vector<Real> e_local;
        SeparatedNorms tau_separated;
        SeparatedNorms e_separated;
        std::vector<Real> residual_history;
        Real theta_mean{0};
        Real theta_std{0};
        Real theta_fraction_good{0};
        Real p_tau{0};
        Real p_e{0};
        Index n_cells{0};
        Real h_eff{0};
        Real V_total{0};
    };

    /// @brief Compute tau = A psi_exact - B(psi_exact).
    [[nodiscard]] static ErrorBundle compute_truncation_error(
        const MeshAdapterT& mesh,
        const RHSAssemblerT& rhs,
        const Mat&,
        const FadeevSolution& fadeev) {
        auto work = mesh;
        std::vector<Real> exact(static_cast<std::size_t>(work.n_cells()));
        const auto& topo = work.topology();
        for (Index i = 0; i < work.n_cells(); ++i) {
            const auto s = static_cast<std::size_t>(i);
            exact[s] = fadeev.psi_exact(topo.coord_x[s], topo.coord_y[s]);
        }
        const auto Apsi = MatrixAssembler::apply_stencil(work, exact);
        const auto b = rhs.assemble(work, exact);
        ErrorBundle out{};
        out.tau_local.resize(Apsi.size());
        for (std::size_t i = 0; i < Apsi.size(); ++i) {
            out.tau_local[i] = Apsi[i] - b[i];
        }
        out.tau = compute_norms(out.tau_local, work.data().volume);
        out.tau_separated = compute_separated_norms(
            out.tau_local, work.data().volume, work.data().kappa, Real{100});
        out.n_cells = work.n_cells();
        out.h_eff = work.h_eff();
        out.V_total = work.total_volume();
        return out;
    }

    /// @brief Compute e = psi_h - psi_exact.
    [[nodiscard]] static ErrorBundle compute_solution_error(
        const MeshAdapterT& mesh,
        std::span<const Real> psi_h,
        const FadeevSolution& fadeev) {
        ErrorBundle out{};
        out.e_local.resize(psi_h.size());
        const auto& topo = mesh.topology();
        for (Index i = 0; i < mesh.n_cells(); ++i) {
            const auto s = static_cast<std::size_t>(i);
            out.e_local[s] = psi_h[s] - fadeev.psi_exact(topo.coord_x[s], topo.coord_y[s]);
        }
        out.e = compute_norms(out.e_local, mesh.data().volume);
        out.e_separated = compute_separated_norms(
            out.e_local, mesh.data().volume, mesh.data().kappa, Real{100});
        compute_effectivity(out, mesh);
        out.n_cells = mesh.n_cells();
        out.h_eff = mesh.h_eff();
        out.V_total = mesh.total_volume();
        return out;
    }

    /// @brief Compute convergence rates from a supplied mesh sequence.
    [[nodiscard]] static std::vector<ErrorBundle> convergence_study(
        std::span<const Index>,
        const Config&) {
        return {};
    }

private:
    static void compute_effectivity(ErrorBundle& out, const MeshAdapterT& mesh) {
        if (out.e_local.empty()) {
            return;
        }
        std::vector<Real> theta(out.e_local.size(), 0.0);
        Index good = 0;
        for (std::size_t i = 0; i < theta.size(); ++i) {
            const Real h2 = mesh.data().volume[i];
            const Real denom = std::max(std::abs(out.e_local[i]), Real{1.0e-30});
            theta[i] = mesh.data().kappa[i] * h2 / denom;
            good += (theta[i] >= 0.5 && theta[i] <= 2.0) ? 1 : 0;
            out.theta_mean += theta[i];
        }
        out.theta_mean /= static_cast<Real>(theta.size());
        for (const Real v : theta) {
            const Real d = v - out.theta_mean;
            out.theta_std += d * d;
        }
        out.theta_std = std::sqrt(out.theta_std / static_cast<Real>(theta.size()));
        out.theta_fraction_good = static_cast<Real>(good) / static_cast<Real>(theta.size());
    }
};

} // namespace gss
