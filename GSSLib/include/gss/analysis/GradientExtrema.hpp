#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

#include "gss/benchmark/FadeevSolution.hpp"
#include "gss/core/MeshAdapter.hpp"

namespace gss {

/// @brief Classification of critical points.
enum class CriticalPointType : std::uint8_t {
    OPoint_max = 0,
    OPoint_min = 1,
    XPoint = 2,
    Unknown = 3
};

/// @brief Structure-of-arrays bundle for detected extrema.
struct ExtremaBundle {
    std::vector<Index> cell_ids;
    std::vector<Real> x;
    std::vector<Real> y;
    std::vector<Real> grad_psi;
    std::vector<Real> grad_F;
    std::vector<Real> source;
    std::vector<Real> kappa;
    std::vector<Real> psi_value;
    std::vector<Real> psi_exact;
    std::vector<Real> local_error;
    std::vector<Real> tau_local;
    std::vector<Real> det_H;
    std::vector<Real> trace_H;
    std::vector<CriticalPointType> point_type;
    Index n_points{0};
};

/// @brief Detects gradient extrema and critical points.
template<typename MeshAdapterT>
struct GradientExtremaDetector {
    /// @brief Find top-k cells by gradient magnitude with spatial filtering.
    [[nodiscard]] static ExtremaBundle find_top_k(
        const MeshAdapterT& mesh,
        const CellData& data,
        Index top_k,
        Real min_distance_h_eff,
        const FadeevSolution* fadeev = nullptr) {
        std::vector<Index> ids(static_cast<std::size_t>(data.n_cells));
        std::iota(ids.begin(), ids.end(), 0);
        std::ranges::sort(ids, [&](Index a, Index b) {
            return grad_norm(data, a) > grad_norm(data, b);
        });
        ExtremaBundle out{};
        const Real min_d = min_distance_h_eff * mesh.h_eff();
        for (const Index id : ids) {
            if (out.n_points >= top_k) {
                break;
            }
            if (is_far_enough(mesh, out, id, min_d)) {
                append(mesh, data, out, id, fadeev);
            }
        }
        return out;
    }

    /// @brief Find cells where gradient magnitude is below threshold.
    [[nodiscard]] static ExtremaBundle find_critical_points(
        const MeshAdapterT& mesh,
        const CellData& data,
        Real grad_threshold,
        const FadeevSolution* fadeev = nullptr) {
        ExtremaBundle out{};
        for (Index i = 0; i < data.n_cells; ++i) {
            if (grad_norm(data, i) < grad_threshold) {
                append(mesh, data, out, i, fadeev);
                const Real det = 0.0;
                const Real trace = 0.0;
                out.det_H.back() = det;
                out.trace_H.back() = trace;
                out.point_type.back() = CriticalPointType::Unknown;
            }
        }
        return out;
    }

    /// @brief Topological count summary for critical points.
    struct TopologyReport {
        Index n_xpoints{0};
        Index n_opoints{0};
        Real mean_grad_at_xpoints{0};
        Real mean_grad_at_opoints{0};
    };

    /// @brief Count X and O points from a critical-point bundle.
    [[nodiscard]] static TopologyReport topology_report(const ExtremaBundle& critical) {
        TopologyReport out{};
        for (std::size_t i = 0; i < critical.point_type.size(); ++i) {
            if (critical.point_type[i] == CriticalPointType::XPoint) {
                ++out.n_xpoints;
                out.mean_grad_at_xpoints += critical.grad_psi[i];
            } else if (critical.point_type[i] == CriticalPointType::OPoint_max
                       || critical.point_type[i] == CriticalPointType::OPoint_min) {
                ++out.n_opoints;
                out.mean_grad_at_opoints += critical.grad_psi[i];
            }
        }
        if (out.n_xpoints > 0) {
            out.mean_grad_at_xpoints /= static_cast<Real>(out.n_xpoints);
        }
        if (out.n_opoints > 0) {
            out.mean_grad_at_opoints /= static_cast<Real>(out.n_opoints);
        }
        return out;
    }

private:
    [[nodiscard]] static Real grad_norm(const CellData& data, Index i) noexcept {
        const auto s = static_cast<std::size_t>(i);
        return std::hypot(data.grad_psi_x[s], data.grad_psi_y[s]);
    }

    [[nodiscard]] static bool is_far_enough(const MeshAdapterT& mesh, const ExtremaBundle& out, Index id, Real min_d) {
        const auto& topo = mesh.topology();
        const auto s = static_cast<std::size_t>(id);
        for (std::size_t i = 0; i < out.x.size(); ++i) {
            if (std::hypot(out.x[i] - topo.coord_x[s], out.y[i] - topo.coord_y[s]) < min_d) {
                return false;
            }
        }
        return true;
    }

    static void append(const MeshAdapterT& mesh, const CellData& data, ExtremaBundle& out, Index id, const FadeevSolution* fadeev) {
        const auto& topo = mesh.topology();
        const auto s = static_cast<std::size_t>(id);
        out.cell_ids.push_back(id);
        out.x.push_back(topo.coord_x[s]);
        out.y.push_back(topo.coord_y[s]);
        out.grad_psi.push_back(grad_norm(data, id));
        out.grad_F.push_back(std::hypot(data.grad_F_x[s], data.grad_F_y[s]));
        out.source.push_back(std::abs(data.source[s]));
        out.kappa.push_back(data.kappa[s]);
        out.psi_value.push_back(data.psi[s]);
        const Real exact = fadeev ? fadeev->psi_exact(topo.coord_x[s], topo.coord_y[s]) : 0.0;
        out.psi_exact.push_back(exact);
        out.local_error.push_back(std::abs(data.psi[s] - exact));
        out.tau_local.push_back(0.0);
        out.det_H.push_back(0.0);
        out.trace_H.push_back(0.0);
        out.point_type.push_back(CriticalPointType::Unknown);
        out.n_points = static_cast<Index>(out.cell_ids.size());
    }
};

} // namespace gss
