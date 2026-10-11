#include <gtest/gtest.h>

#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/core/CellPermutation.hpp"
#include "gss/linear/PETScLinearSystem.hpp"

TEST(MatrixAssembler, ConstantFieldHasZeroStencil) {
    gss::Config config{};
    config.mesh.n_generators = 64;
    const auto mesh = gss::MeshAdapter::uniform_from_config(config);
    std::vector<gss::Real> one(static_cast<std::size_t>(mesh.n_cells()), 1.0);
    const auto applied = gss::MatrixAssembler::apply_stencil(mesh, one);
    for (const auto v : applied) {
        EXPECT_NEAR(v, 0.0, 1.0e-12);
    }
}

TEST(MatrixAssembler, MatrixFreeStencilIsEquivalentAfterPermutation) {
    gss::Config config{};
    config.mesh.n_generators = 16;
    const auto mesh = gss::MeshAdapter::uniform_from_config(config);

    std::vector<gss::Real> psi(static_cast<std::size_t>(mesh.n_cells()), 0.0);
    for (gss::Index i = 0; i < mesh.n_cells(); ++i) {
        const auto s = static_cast<std::size_t>(i);
        psi[s] = mesh.topology().coord_x[s] + 2.0 * mesh.topology().coord_y[s];
    }
    const auto reference = gss::MatrixAssembler::apply_stencil(mesh, psi);

    auto topology = mesh.topology();
    auto volumes = mesh.data().volume;
    const auto permutation = gss::CellPermutation::from_new_to_old(
        {5, 0, 7, 2, 10, 1, 12, 3, 14, 4, 15, 6, 13, 8, 11, 9});
    gss::apply_cell_permutation(topology, volumes, permutation);
    const gss::MeshAdapter permuted_mesh{std::move(topology), std::move(volumes)};

    std::vector<gss::Real> permuted_psi(static_cast<std::size_t>(permuted_mesh.n_cells()), 0.0);
    for (gss::Index new_index = 0; new_index < permuted_mesh.n_cells(); ++new_index) {
        const auto ns = static_cast<std::size_t>(new_index);
        const auto old_index = permutation.new_to_old()[ns];
        permuted_psi[ns] = psi[static_cast<std::size_t>(old_index)];
    }

    const auto permuted_result = gss::MatrixAssembler::apply_stencil(permuted_mesh, permuted_psi);
    for (gss::Index new_index = 0; new_index < permuted_mesh.n_cells(); ++new_index) {
        const auto ns = static_cast<std::size_t>(new_index);
        const auto old_index = permutation.new_to_old()[ns];
        EXPECT_NEAR(
            permuted_result[ns],
            reference[static_cast<std::size_t>(old_index)],
            1.0e-12);
    }
}
