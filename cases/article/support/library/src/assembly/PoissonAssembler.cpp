#include "gsmaker/assembly/PoissonAssembler.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>





namespace gss {
namespace {

[[nodiscard]] Real denominator(const EvaluatedBoundaryFace& face) {
    const auto& c = face.coefficients;
    const auto& g = face.geometry;
    const Real value = c.alpha * g.perpendicular_distance + c.beta;
    const Real scale = std::abs(c.alpha) * g.perpendicular_distance + std::abs(c.beta);
    const Real tolerance = Real{64} * std::numeric_limits<Real>::epsilon() * std::max(scale, Real{1});
    if (!(std::abs(value) > tolerance)) {
        throw std::runtime_error("singular boundary equation: alpha*d + beta is zero");
    }
    return value;
}

void validate_diffusivity(Real diffusivity) {
    if (!std::isfinite(static_cast<double>(diffusivity)) || !(diffusivity > Real{0})) {
        throw std::invalid_argument("Poisson diffusivity must be finite and positive");
    }
}

} // namespace

PoissonSystem PoissonAssembler::assemble(
    const MeshAdapter& mesh,
    std::span<const EvaluatedBoundaryFace> boundary,
    Real diffusivity) {

    validate_diffusivity(diffusivity);
    const Index n = mesh.n_cells();
    const auto ns = static_cast<std::size_t>(n);
    const auto& topo = mesh.topology();

    std::vector<std::map<Index, Real>> rows(ns);
    std::vector<Real> boundary_rhs(ns, Real{0});

    for (Index p = 0; p < n; ++p) {
        const auto ps = static_cast<std::size_t>(p);
        Real diagonal = Real{0};
        for (Index e = topo.neighbor_start[ps]; e < topo.neighbor_start[ps + 1U]; ++e) {
            const auto es = static_cast<std::size_t>(e);
            const Index q = topo.neighbor_list[es];
            const Real coefficient = diffusivity * topo.interface_area[es] / topo.distance[es];
            diagonal += coefficient;
            rows[ps][q] -= coefficient;
        }
        rows[ps][p] += diagonal;
    }

    for (const auto& face : boundary) {
        const Index p = face.geometry.owner_cell;
        if (p < 0 || p >= n) throw std::runtime_error("boundary contribution owner is out of range");
        const Real d = denominator(face);
        const Real measure = face.geometry.measure;
        const Real diagonal = diffusivity * measure * face.coefficients.alpha / d;
        const Real rhs = diffusivity * measure * face.coefficients.gamma / d;
        rows[static_cast<std::size_t>(p)][p] += diagonal;
        boundary_rhs[static_cast<std::size_t>(p)] += rhs;
    }

    std::vector<Index> row_start;
    std::vector<Index> column;
    std::vector<Real> value;
    row_start.reserve(ns + 1U);
    row_start.push_back(Index{0});

    for (Index p = 0; p < n; ++p) {
        for (const auto& [q, a] : rows[static_cast<std::size_t>(p)]) {
            if (a == Real{0}) continue;
            column.push_back(q);
            value.push_back(a);
        }
        row_start.push_back(static_cast<Index>(column.size()));
    }

    CsrMatrixProperties properties;

    // The integrated two-point diffusion operator is symmetric by
    // construction.  Each internal face uses the reciprocal transmissibility
    // T_pq = |f_pq| / d_pq = T_qp; boundary conditions modify only the owner
    // diagonal.  Positive definiteness is deliberately not asserted here,
    // because it depends on the boundary-value problem (for example, pure
    // homogeneous Neumann data leave a constant nullspace).
    properties.symmetric = true;
    properties.structurally_symmetric = true;

    return PoissonSystem{
        CsrMatrix{
            n,
            std::move(row_start),
            std::move(column),
            std::move(value),
            properties},
        std::move(boundary_rhs)};
}

std::vector<Real> PoissonAssembler::build_rhs(
    const MeshAdapter& mesh,
    std::span<const Real> source,
    std::span<const EvaluatedBoundaryFace> boundary,
    Real diffusivity) {

    if (source.size() != static_cast<std::size_t>(mesh.n_cells())) {
        throw std::invalid_argument("Poisson source array size mismatch");
    }
    const auto system = assemble(mesh, boundary, diffusivity);
    std::vector<Real> rhs = system.boundary_rhs;
    const auto volumes = mesh.geometry().cell_volume();
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        if (!std::isfinite(static_cast<double>(source[i]))) {
            throw std::invalid_argument("Poisson source contains a non-finite value");
        }
        rhs[i] -= volumes[i] * source[i];
    }
    return rhs;
}

} // namespace gss
