#include <gtest/gtest.h>

#include "gss/analysis/StructuralMetrics.hpp"
#include "gss/core/MeshAdapter.hpp"

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

} // namespace

TEST(StructuralMetrics, ComputesChainMetrics) {
    const auto metrics = gss::compute_structural_metrics(make_three_cell_chain());

    EXPECT_EQ(metrics.bandwidth, 1);
    EXPECT_EQ(metrics.profile, 2);
    EXPECT_EQ(metrics.edge_count, 2);
    EXPECT_EQ(metrics.matrix_nnz, 7);
    EXPECT_EQ(metrics.lower_triangle_nnz, 5);
    EXPECT_DOUBLE_EQ(metrics.adjacency_index_span.mean, 1.0);
    EXPECT_DOUBLE_EQ(metrics.adjacency_index_span.median, 1.0);
    EXPECT_DOUBLE_EQ(metrics.adjacency_index_span.maximum, 1.0);
    EXPECT_DOUBLE_EQ(metrics.adjacency_index_span.standard_deviation, 0.0);
}

TEST(StructuralMetrics, DetectsUniformGridBandwidth) {
    gss::Config config{};
    config.mesh.n_generators = 16;
    const auto mesh = gss::MeshAdapter::uniform_from_config(config);

    const auto metrics = gss::compute_structural_metrics(mesh.topology());

    EXPECT_EQ(metrics.bandwidth, 4);
    EXPECT_EQ(metrics.edge_count, 24);
    EXPECT_EQ(metrics.matrix_nnz, 64);
    EXPECT_EQ(metrics.lower_triangle_nnz, 40);
    EXPECT_GT(metrics.profile, 0);
    EXPECT_GT(metrics.adjacency_index_span.mean, 0.0);
    EXPECT_GE(metrics.adjacency_index_span.maximum, metrics.adjacency_index_span.mean);
}
