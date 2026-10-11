#include <gtest/gtest.h>

#include "gss/analysis/GradientExtrema.hpp"

TEST(GradientExtrema, BundleUsesSeparateArrays) {
    gss::Config config{};
    config.mesh.n_generators = 64;
    auto mesh = gss::MeshAdapter::uniform_from_config(config);
    mesh.data().grad_psi_x.assign(static_cast<std::size_t>(mesh.n_cells()), 1.0);
    mesh.data().grad_psi_y.assign(static_cast<std::size_t>(mesh.n_cells()), 0.0);
    const auto out = gss::GradientExtremaDetector<gss::MeshAdapter>::find_top_k(mesh, mesh.data(), 5, 1.0);
    EXPECT_LE(out.n_points, 5);
    EXPECT_EQ(out.cell_ids.size(), out.x.size());
    EXPECT_EQ(out.x.size(), out.y.size());
}
