#include <gtest/gtest.h>

#include "gss/assembly/RHSAssembler.hpp"

TEST(RHSAssembler, ProducesVolumeWeightedSource) {
    gss::Config config{};
    config.mesh.n_generators = 64;
    auto mesh = gss::MeshAdapter::uniform_from_config(config);
    std::vector<gss::Real> psi(static_cast<std::size_t>(mesh.n_cells()), 0.0);
    gss::RHSAssembler<> rhs(config);
    const auto b = rhs.assemble(mesh, psi);
    ASSERT_EQ(b.size(), psi.size());
    for (std::size_t i = 0; i < b.size(); ++i) {
        EXPECT_DOUBLE_EQ(b[i], mesh.data().source[i] * mesh.data().volume[i]);
    }
}
