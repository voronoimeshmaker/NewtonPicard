#pragma once

#include <algorithm>
#include <cmath>
#include <span>
#include <utility>
#include <vector>

#include "gss/benchmark/FadeevSolution.hpp"
#include "gss/core/Config.hpp"
#include "gss/core/MeshAdapter.hpp"
#include "gss/linear/PETScLinearSystem.hpp"

namespace gss {

/// @brief Default source policy used for standalone benchmarks.
struct FadeevSourcePolicy {
    FadeevSolution fadeev;

    /// @brief Compute the exact source at a cell center.
    template<typename Scalar>
    [[nodiscard]] Scalar compute_source(Index, Real x, Real y, std::span<const Scalar>) const {
        return fadeev.source_exact(x, y);
    }
};

/// @brief Assemble nonlinear RHS B(psi) = S(psi) * volume.
template<typename SourcePolicy = FadeevSourcePolicy>
class RHSAssembler {
public:
    /// @brief Construct from configuration and source policy.
    explicit RHSAssembler(const Config& config, SourcePolicy source_policy = {})
        : config_(config), source_policy_(std::move(source_policy)) {}

    /// @brief Pure scalar-generic RHS entry, shared by evaluation and ADOL-C recording.
    template<typename Scalar>
    [[nodiscard]] Scalar value_at(const MeshAdapter& mesh, std::span<const Scalar> psi, Index p) const {
        const auto s = static_cast<std::size_t>(p);
        const auto& topo = mesh.topology();
        return source_policy_.compute_source(p, topo.coord_x[s], topo.coord_y[s], psi)
            * mesh.data().volume[s];
    }

    /// @brief Assemble the source and right-hand side arrays.
    [[nodiscard]] std::vector<Real> assemble(MeshAdapter& mesh, std::span<const Real> psi) const {
        auto& data = mesh.data();
        const auto& topo = mesh.topology();
        std::vector<Real> b(static_cast<std::size_t>(data.n_cells), 0.0);
        Real* const source = data.source.data();
        Real* const grad_x = data.grad_psi_x.data();
        Real* const grad_y = data.grad_psi_y.data();
        Real* const kappa = data.kappa.data();
        const Real* const volume = data.volume.data();
        const Real* const x = topo.coord_x.data();
        const Real* const y = topo.coord_y.data();
        const Index n = data.n_cells;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(n > 1024)
#endif
        for (Index p = 0; p < n; ++p) {
            source[p] = source_policy_.compute_source(p, x[p], y[p], psi);
            b[static_cast<std::size_t>(p)] = source[p] * volume[p];
            estimate_gradient_component(mesh, psi, p, grad_x[p], grad_y[p]);
            const Real g2 = grad_x[p] * grad_x[p] + grad_y[p] * grad_y[p];
            kappa[p] = (g2 < config_.gradient_analysis.grad_critical_threshold
                         * config_.gradient_analysis.grad_critical_threshold)
                ? config_.benchmark.kappa_island_threshold * 10.0
                : 1.0;
        }
        return b;
    }

    /// @brief Assemble B(psi) directly into an already allocated distributed PETSc Vec.
    void assemble_into(MeshAdapter& mesh, Vec psi, Vec b) const {
        (void)psi;
        auto& data = mesh.data();
        const auto& topo = mesh.topology();
        PetscInt rstart = 0;
        PetscInt rend = 0;
        petsc_check(VecGetOwnershipRange(b, &rstart, &rend), "get RHS ownership range");
        const Real* const volume = data.volume.data();
        const Real* const x = topo.coord_x.data();
        const Real* const y = topo.coord_y.data();
        Real* const source = data.source.data();
        for (PetscInt i = rstart; i < rend; ++i) {
            const auto s = static_cast<std::size_t>(i);
            source[s] = source_policy_.compute_source(i, x[s], y[s], std::span<const Real>{});
            const Real value = source[s] * volume[s];
            petsc_check(VecSetValue(b, i, value, INSERT_VALUES), "set RHS value");
        }
        petsc_check(VecAssemblyBegin(b), "begin RHS assembly");
        petsc_check(VecAssemblyEnd(b), "end RHS assembly");
    }

private:
    static void estimate_gradient_component(
        const MeshAdapter& mesh,
        std::span<const Real> psi,
        Index p,
        Real& gx,
        Real& gy) noexcept {
        const auto& topo = mesh.topology();
        const Index begin = topo.neighbor_start[static_cast<std::size_t>(p)];
        const Index end = topo.neighbor_start[static_cast<std::size_t>(p + 1)];
        gx = 0.0;
        gy = 0.0;
        Real wsum = 0.0;
        for (Index e = begin; e < end; ++e) {
            const Index q = topo.neighbor_list[static_cast<std::size_t>(e)];
            const Real dx = topo.coord_x[static_cast<std::size_t>(q)] - topo.coord_x[static_cast<std::size_t>(p)];
            const Real dy = topo.coord_y[static_cast<std::size_t>(q)] - topo.coord_y[static_cast<std::size_t>(p)];
            const Real d2 = std::max(dx * dx + dy * dy, Real{1.0e-30});
            const Real diff = psi[static_cast<std::size_t>(q)] - psi[static_cast<std::size_t>(p)];
            gx += diff * dx / d2;
            gy += diff * dy / d2;
            wsum += 1.0;
        }
        if (wsum > 0.0) {
            gx /= wsum;
            gy /= wsum;
        }
    }

    Config config_;
    SourcePolicy source_policy_;
};

} // namespace gss
