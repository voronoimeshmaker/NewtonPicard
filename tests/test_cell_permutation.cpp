#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "gss/core/CellPermutation.hpp"

namespace {

[[nodiscard]] gss::MeshTopology make_three_cell_chain() {
    gss::MeshTopology topology;
    topology.coord_x = {0.0, 1.0, 2.0};
    topology.coord_y = {0.0, 0.0, 0.0};
    topology.is_boundary_cell = {1U, 0U, 1U};
    topology.neighbor_start = {0, 1, 3, 4};
    topology.neighbor_list = {1, 0, 2, 1};
    topology.interface_area = {1.0, 1.0, 1.0, 1.0};
    topology.distance = {1.0, 1.0, 1.0, 1.0};
    return topology;
}

[[nodiscard]] std::vector<std::pair<gss::Index, gss::Index>>
undirected_edges(const gss::MeshTopology& topology) {
    std::vector<std::pair<gss::Index, gss::Index>> edges;
    for (gss::Index p = 0; p + 1 < static_cast<gss::Index>(topology.neighbor_start.size()); ++p) {
        const gss::Index begin = topology.neighbor_start[static_cast<std::size_t>(p)];
        const gss::Index end = topology.neighbor_start[static_cast<std::size_t>(p + 1)];
        for (gss::Index edge = begin; edge < end; ++edge) {
            const gss::Index q = topology.neighbor_list[static_cast<std::size_t>(edge)];
            edges.emplace_back(std::min(p, q), std::max(p, q));
        }
    }
    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

} // namespace

TEST(CellPermutation, RejectsNonBijectiveMappings) {
    EXPECT_THROW(
        static_cast<void>(gss::CellPermutation::from_new_to_old({0, 0, 2})),
        std::runtime_error);
    EXPECT_THROW(
        static_cast<void>(gss::CellPermutation::from_new_to_old({0, 1, 3})),
        std::runtime_error);
}

TEST(CellPermutation, BuildsInverseAndComposition) {
    const auto permutation = gss::CellPermutation::from_new_to_old({2, 0, 1});
    const auto inverse = permutation.inverse();

    ASSERT_EQ(inverse.new_to_old().size(), 3U);
    EXPECT_EQ(inverse.new_to_old()[0], 1);
    EXPECT_EQ(inverse.new_to_old()[1], 2);
    EXPECT_EQ(inverse.new_to_old()[2], 0);

    const auto identity = permutation.then(inverse);
    EXPECT_TRUE(identity.is_identity());
}

TEST(CellPermutation, AppliesToTopologyAndVolumes) {
    auto topology = make_three_cell_chain();
    std::vector<gss::Real> volumes{10.0, 20.0, 30.0};
    const auto original_edges = undirected_edges(topology);

    const auto permutation = gss::CellPermutation::from_new_to_old({2, 0, 1});
    gss::apply_cell_permutation(topology, volumes, permutation);

    EXPECT_EQ(topology.coord_x[0], 2.0);
    EXPECT_EQ(topology.coord_x[1], 0.0);
    EXPECT_EQ(topology.coord_x[2], 1.0);
    EXPECT_EQ(volumes[0], 30.0);
    EXPECT_EQ(volumes[1], 10.0);
    EXPECT_EQ(volumes[2], 20.0);
    EXPECT_EQ(topology.is_boundary_cell[0], 1U);
    EXPECT_EQ(topology.is_boundary_cell[1], 1U);
    EXPECT_EQ(topology.is_boundary_cell[2], 0U);

    std::vector<std::pair<gss::Index, gss::Index>> restored_edges;
    for (const auto [p, q] : undirected_edges(topology)) {
        const auto old_p = permutation.new_to_old()[static_cast<std::size_t>(p)];
        const auto old_q = permutation.new_to_old()[static_cast<std::size_t>(q)];
        restored_edges.emplace_back(std::min(old_p, old_q), std::max(old_p, old_q));
    }
    std::sort(restored_edges.begin(), restored_edges.end());
    EXPECT_EQ(restored_edges, original_edges);
}
