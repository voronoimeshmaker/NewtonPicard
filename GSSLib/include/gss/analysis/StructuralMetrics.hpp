#pragma once

#include "gss/core/Types.hpp"

namespace gss {

/// @brief Summary statistics for adjacency index spans.
struct AdjacencySpanStats {
    Real mean{0};
    Real median{0};
    Real p90{0};
    Real p95{0};
    Real p99{0};
    Real maximum{0};
    Real standard_deviation{0};
};

/// @brief Structural metrics of the symmetric cell adjacency graph.
struct StructuralMetrics {
    Index bandwidth{0};
    long long profile{0};
    Index edge_count{0};
    long long matrix_nnz{0};
    long long lower_triangle_nnz{0};
    AdjacencySpanStats adjacency_index_span{};
};

/// @brief Compute graph bandwidth, lower profile, and adjacency-span statistics.
[[nodiscard]] StructuralMetrics compute_structural_metrics(const MeshTopology& topology);

} // namespace gss
