#include <gtest/gtest.h>

#include <limits>

#include "gss/assembly/MatrixAssembler.hpp"

TEST(PatchTest, LinearFieldInteriorResidualIsSmall) {
    gss::Config config{};
    config.mesh.n_generators = 64;
    const auto mesh = gss::MeshAdapter::uniform_from_config(config);
    std::vector<gss::Real> psi(static_cast<std::size_t>(mesh.n_cells()));
    for (gss::Index i = 0; i < mesh.n_cells(); ++i) {
        const auto s = static_cast<std::size_t>(i);
        psi[s] = 2.0 * mesh.topology().coord_x[s] + 3.0 * mesh.topology().coord_y[s] + 1.0;
    }
    const auto tau = gss::MatrixAssembler::apply_stencil(mesh, psi);
    gss::Real max_abs = 0.0;
    for (const auto v : tau) {
        max_abs = std::max(max_abs, std::abs(v));
    }
    EXPECT_LT(max_abs, 20.0);
}
