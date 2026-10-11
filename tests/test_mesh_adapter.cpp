#include <gtest/gtest.h>

#include "gss/core/MeshAdapter.hpp"

TEST(MeshAdapter, UniformMeshHasConsistentCSRAndPositiveMetrics) {
    gss::Config config{};
    config.mesh.n_generators = 64;
    const auto mesh = gss::MeshAdapter::uniform_from_config(config);
    EXPECT_EQ(mesh.n_cells(), 64);
    EXPECT_EQ(mesh.topology().neighbor_start.back(), static_cast<gss::Index>(mesh.topology().neighbor_list.size()));
    for (const auto v : mesh.data().volume) {
        EXPECT_GT(v, 0.0);
    }
    for (const auto d : mesh.topology().distance) {
        EXPECT_GT(d, 0.0);
    }
    for (const auto a : mesh.topology().interface_area) {
        EXPECT_GT(a, 0.0);
    }
}
