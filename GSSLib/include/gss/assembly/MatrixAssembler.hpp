#pragma once

#include <stdexcept>
#include <vector>

#include <petsc.h>

#include "gss/core/MeshAdapter.hpp"

namespace gss {

/// @brief Assembles the FVM Voronoi diffusion operator into a PETSc matrix.
class MatrixAssembler {
public:
    template<typename Scalar>
    [[nodiscard]] static Scalar stencil_at(const MeshAdapter& mesh, std::span<const Scalar> psi,
                                           Index p, Real diffusivity = Real{1}) {
        const auto& topo = mesh.topology();
        Scalar sum = 0.0;
        const auto ps = static_cast<std::size_t>(p);
        for (Index e = topo.neighbor_start[ps]; e < topo.neighbor_start[ps + 1]; ++e) {
            const auto es = static_cast<std::size_t>(e);
            const auto q = static_cast<std::size_t>(topo.neighbor_list[es]);
            const Real c = diffusivity * topo.interface_area[es] / topo.distance[es];
            sum += c * (psi[ps] - psi[q]);
        }
        return sum;
    }

    /// @brief Assemble an SPD diffusion matrix from mesh topology.
    [[nodiscard]] static Mat assemble(const MeshAdapter& mesh, MPI_Comm comm = PETSC_COMM_WORLD,
                                      Real diffusivity = Real{1});

    /// @brief Apply the matrix-free FVM stencil to a PETSc vector.
    static void apply_stencil(const MeshAdapter& mesh, Vec psi, Vec out,
                              Real diffusivity = Real{1});

    /// @brief Apply the matrix-free FVM stencil to a vector.
    [[nodiscard]] static std::vector<Real> apply_stencil(const MeshAdapter& mesh, std::span<const Real> psi,
                                                         Real diffusivity = Real{1});
};

} // namespace gss
