#include "gss/analysis/StructuralMetrics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace gss {
namespace {

[[nodiscard]] Real percentile(std::vector<Real>& values, Real fraction) {
    if (values.empty()) {
        return Real{0};
    }
    fraction = std::clamp(fraction, Real{0}, Real{1});
    const Real scaled = fraction * static_cast<Real>(values.size() - 1U);
    const auto index = static_cast<std::size_t>(std::llround(scaled));
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(index), values.end());
    return values[index];
}

} // namespace

StructuralMetrics compute_structural_metrics(const MeshTopology& topology) {
    const auto n = static_cast<Index>(topology.coord_x.size());
    if (topology.neighbor_start.size() != static_cast<std::size_t>(n + 1)) {
        throw std::runtime_error("cannot compute structural metrics for inconsistent topology");
    }

    StructuralMetrics metrics{};
    std::vector<Real> spans;
    spans.reserve(topology.neighbor_list.size() / 2U);

    for (Index row = 0; row < n; ++row) {
        const Index begin = topology.neighbor_start[static_cast<std::size_t>(row)];
        const Index end = topology.neighbor_start[static_cast<std::size_t>(row + 1)];
        Index min_column = row;
        for (Index edge = begin; edge < end; ++edge) {
            const Index column = topology.neighbor_list[static_cast<std::size_t>(edge)];
            if (column < 0 || column >= n) {
                throw std::runtime_error("structural metrics found an out-of-range neighbour");
            }
            const Index span = std::abs(row - column);
            metrics.bandwidth = std::max(metrics.bandwidth, span);
            min_column = std::min(min_column, column);
            if (row < column) {
                spans.push_back(static_cast<Real>(span));
            }
        }
        metrics.profile += static_cast<long long>(row - min_column);
    }

    metrics.edge_count = static_cast<Index>(spans.size());
    metrics.matrix_nnz = static_cast<long long>(n) + 2LL * static_cast<long long>(metrics.edge_count);
    metrics.lower_triangle_nnz = static_cast<long long>(n) + static_cast<long long>(metrics.edge_count);
    if (spans.empty()) {
        return metrics;
    }

    const Real sum = std::accumulate(spans.begin(), spans.end(), Real{0});
    metrics.adjacency_index_span.mean = sum / static_cast<Real>(spans.size());
    metrics.adjacency_index_span.maximum = *std::max_element(spans.begin(), spans.end());

    Real variance = 0;
    for (const Real value : spans) {
        const Real delta = value - metrics.adjacency_index_span.mean;
        variance += delta * delta;
    }
    variance /= static_cast<Real>(spans.size());
    metrics.adjacency_index_span.standard_deviation = std::sqrt(variance);

    auto work = spans;
    metrics.adjacency_index_span.median = percentile(work, Real{0.50});
    work = spans;
    metrics.adjacency_index_span.p90 = percentile(work, Real{0.90});
    work = spans;
    metrics.adjacency_index_span.p95 = percentile(work, Real{0.95});
    work = spans;
    metrics.adjacency_index_span.p99 = percentile(work, Real{0.99});
    return metrics;
}

} // namespace gss
