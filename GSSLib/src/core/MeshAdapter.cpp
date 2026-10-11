#include "gss/core/MeshAdapter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <utility>

#include "gss/core/CellPermutation.hpp"
#include "gss/core/Reordering.hpp"

#include <VoronoiMeshMaker/Boundary2D/Builders/Boundary2DBuilder.hpp>
#include <VoronoiMeshMaker/Boundary2D/Shapes/Rectangle.hpp>
#include <VoronoiMeshMaker/Sites2D/Factory/SiteFactory.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiBuilder2D.hpp>

namespace gss {
namespace {

[[nodiscard]] Real distance_between(Real ax, Real ay, Real bx, Real by) noexcept {
    return std::hypot(bx - ax, by - ay);
}

template<typename Diagram, typename Cell, typename Edge>
[[nodiscard]] Index neighbor_volume_for_edge(
    const Diagram& diagram,
    const Cell& cell,
    const Edge& edge,
    Real px,
    Real py) {
    constexpr Real tol_scale = Real{1.0e-8};
    const Real center_dist2 =
        (edge.midpoint.x - px) * (edge.midpoint.x - px) +
        (edge.midpoint.y - py) * (edge.midpoint.y - py);

    Real best_residual = std::numeric_limits<Real>::max();
    Index best_volume = -1;
    for (const auto neighbor_site_id : cell.neighbor_ids) {
        const auto neighbor_site_index = static_cast<std::size_t>(neighbor_site_id.value);
        const auto& neighbor_site = diagram.sites[neighbor_site_index];
        const Real neighbor_dist2 =
            (edge.midpoint.x - neighbor_site.point.x) * (edge.midpoint.x - neighbor_site.point.x) +
            (edge.midpoint.y - neighbor_site.point.y) * (edge.midpoint.y - neighbor_site.point.y);
        const Real residual = std::abs(center_dist2 - neighbor_dist2);
        if (residual < best_residual) {
            best_residual = residual;
            best_volume = static_cast<Index>(diagram.cell(neighbor_site_id).volume_id);
        }
    }

    const Real tol = tol_scale * std::max(Real{1}, center_dist2);
    return best_residual <= tol ? best_volume : Index{-1};
}

/// @brief Compute the cell renumbering permutation for a mesh: perm[new] = old.
[[nodiscard]] std::vector<Index> make_ordering(const Config& config, const MeshTopology& topo) {
    const auto n = static_cast<std::size_t>(topo.coord_x.size());
    std::vector<Index> perm(n);
    std::iota(perm.begin(), perm.end(), Index{0});
    const std::string& r = config.mesh.reordering;

    if (r.empty() || r == "natural") {
        return perm;
    }
    if (r == "random") {
        std::mt19937 rng(static_cast<std::uint32_t>(config.mesh.random_seed));
        std::shuffle(perm.begin(), perm.end(), rng);
        return perm;
    }
    if (r == "morton" || r == "hilbert") {
        const auto xminmax = std::minmax_element(topo.coord_x.begin(), topo.coord_x.end());
        const auto yminmax = std::minmax_element(topo.coord_y.begin(), topo.coord_y.end());
        const Real xmin = *xminmax.first;
        const Real xspan = *xminmax.second - xmin;
        const Real ymin = *yminmax.first;
        const Real yspan = *yminmax.second - ymin;
        const std::uint32_t side = std::uint32_t{1} << reorder::kSfcBits;
        std::vector<std::uint64_t> keys(n);
        for (std::size_t i = 0; i < n; ++i) {
            const std::uint32_t xi = reorder::quantize(topo.coord_x[i], xmin, xspan);
            const std::uint32_t yi = reorder::quantize(topo.coord_y[i], ymin, yspan);
            keys[i] = (r == "morton")
                ? reorder::morton_key(xi, yi)
                : reorder::hilbert_key(xi, yi, side);
        }
        std::sort(perm.begin(), perm.end(), [&](Index a, Index b) {
            return keys[static_cast<std::size_t>(a)] < keys[static_cast<std::size_t>(b)];
        });
        return perm;
    }
    if (r == "rcm") {
        // Reverse Cuthill-McKee on the dual graph (cells = vertices, faces = edges).
        const Index nv = static_cast<Index>(topo.coord_x.size());
        std::vector<std::size_t> deg(static_cast<std::size_t>(nv));
        for (Index v = 0; v < nv; ++v) {
            deg[static_cast<std::size_t>(v)] = static_cast<std::size_t>(
                topo.neighbor_start[static_cast<std::size_t>(v + 1)]
                - topo.neighbor_start[static_cast<std::size_t>(v)]);
        }
        // Pseudo-peripheral start: lowest degree, ties -> lowest index.
        std::size_t start = 0;
        for (std::size_t v = 1; v < deg.size(); ++v) {
            if (deg[v] < deg[start]) start = v;
        }
        std::vector<bool> visited(static_cast<std::size_t>(nv), false);
        std::vector<std::size_t> order;
        order.reserve(static_cast<std::size_t>(nv));
        std::vector<std::size_t> level;
        level.reserve(static_cast<std::size_t>(nv));
        visited[start] = true;
        level.push_back(start);
        while (!level.empty()) {
            // Sort the current level by ascending degree (ties by index).
            std::sort(level.begin(), level.end(), [&](std::size_t a, std::size_t b) {
                return deg[a] != deg[b] ? deg[a] < deg[b] : a < b;
            });
            for (const auto v : level) order.push_back(v);
            std::vector<std::size_t> next;
            next.reserve(level.size() * 4);
            for (const auto v : level) {
                const Index b = topo.neighbor_start[v];
                const Index e = topo.neighbor_start[v + 1];
                for (Index ei = b; ei < e; ++ei) {
                    const auto w = static_cast<std::size_t>(topo.neighbor_list[static_cast<std::size_t>(ei)]);
                    if (!visited[w]) {
                        visited[w] = true;
                        next.push_back(w);
                    }
                }
            }
            level.swap(next);
        }
        // Reverse the BFS order -> RCM.
        std::reverse(order.begin(), order.end());
        for (std::size_t newpos = 0; newpos < order.size(); ++newpos) {
            perm[newpos] = static_cast<Index>(order[newpos]);
        }
        return perm;
    }
    throw std::runtime_error("unknown mesh.reordering: " + r);
}

/// @brief Reorder a freshly built mesh according to the config and record the permutation.
[[nodiscard]] std::vector<Index> reorder_mesh(const Config& config, MeshTopology& topo, std::vector<Real>& volumes) {
    std::vector<Index> perm = make_ordering(config, topo);
    if (!(perm.empty())) {
        bool identity = true;
        for (std::size_t i = 0; i < perm.size(); ++i) {
            if (perm[i] != static_cast<Index>(i)) {
                identity = false;
                break;
            }
        }
        if (!identity) {
            apply_cell_permutation(topo, volumes, CellPermutation::from_new_to_old(perm));
        } else {
            perm.clear();
        }
    }
    return perm;
}

} // namespace

MeshAdapter::MeshAdapter(MeshTopology topology, std::vector<Real> volumes)
    : topology_(std::move(topology)) {
    data_.resize(static_cast<Index>(volumes.size()));
    data_.volume = std::move(volumes);
    validate();
}

MeshAdapter MeshAdapter::uniform_from_config(const Config& config) {
    const Index n = config.mesh.n_generators;
    const Index nx = static_cast<Index>(std::sqrt(static_cast<Real>(n)));
    const Index ny = (nx * nx == n) ? nx : static_cast<Index>(std::ceil(static_cast<Real>(n) / static_cast<Real>(nx)));
    const Real xmin = config.mesh.domain[0];
    const Real xmax = config.mesh.domain[1];
    const Real ymin = config.mesh.domain[2];
    const Real ymax = config.mesh.domain[3];
    const Real dx = (xmax - xmin) / static_cast<Real>(nx);
    const Real dy = (ymax - ymin) / static_cast<Real>(ny);

    MeshTopology topo;
    topo.neighbor_start.reserve(static_cast<std::size_t>(n) + 1);
    topo.coord_x.resize(static_cast<std::size_t>(n));
    topo.coord_y.resize(static_cast<std::size_t>(n));
    topo.is_boundary_cell.assign(static_cast<std::size_t>(n), 0U);
    std::vector<Real> volumes(static_cast<std::size_t>(n), dx * dy);
    topo.neighbor_start.push_back(0);
    for (Index j = 0; j < ny; ++j) {
        for (Index i = 0; i < nx; ++i) {
            const Index p = j * nx + i;
            if (p >= n) {
                break;
            }
            const auto ps = static_cast<std::size_t>(p);
            topo.coord_x[ps] = xmin + (static_cast<Real>(i) + 0.5) * dx;
            topo.coord_y[ps] = ymin + (static_cast<Real>(j) + 0.5) * dy;
            topo.is_boundary_cell[ps] = static_cast<unsigned char>(
                i == 0 || j == 0 || i == nx - 1 || j == ny - 1);
            auto add_neighbor = [&](bool present, Index q, Real area, Real dist) {
                if (!present) {
                    return;
                }
                topo.neighbor_list.push_back(q);
                topo.interface_area.push_back(area);
                topo.distance.push_back(dist);
            };
            add_neighbor(i > 0, p - 1, dy, dx);
            add_neighbor(i + 1 < nx && p + 1 < n, p + 1, dy, dx);
            add_neighbor(j > 0, p - nx, dx, dy);
            add_neighbor(j + 1 < ny && p + nx < n, p + nx, dx, dy);
            topo.neighbor_start.push_back(static_cast<Index>(topo.neighbor_list.size()));
        }
    }
    std::vector<Index> perm = reorder_mesh(config, topo, volumes);
    MeshAdapter mesh{std::move(topo), std::move(volumes)};
    mesh.set_ordering(std::move(perm));
    return mesh;
}

MeshAdapter MeshAdapter::from_config(const Config& config) {
    if (config.mesh.type == "uniform_internal") {
        return uniform_from_config(config);
    }

    namespace b2d = vmm::b2d;
    namespace s2d = vmm::s2d;
    namespace vd2d = vmm::vd2d;

    const Real xmin = config.mesh.domain[0];
    const Real xmax = config.mesh.domain[1];
    const Real ymin = config.mesh.domain[2];
    const Real ymax = config.mesh.domain[3];
    const Real width = xmax - xmin;
    const Real height = ymax - ymin;
    if (width <= 0.0 || height <= 0.0) {
        throw std::runtime_error("mesh.domain must have xmax > xmin and ymax > ymin");
    }

    const auto boundary = b2d::make_boundary(
        b2d::Rectangle(b2d::Point2{xmin, ymin}, width, height));

    const auto n_dir = static_cast<std::size_t>(std::max<Index>(config.mesh.n_per_direction, 1));
    s2d::SiteSet sites;
    if (config.mesh.type == "cartesian_square" || config.mesh.type == "voronoi_uniform") {
        const auto nx_cells = static_cast<std::size_t>(config.mesh.n_per_direction);
        const auto ny_cells = config.mesh.n_per_direction_y > 0
            ? static_cast<std::size_t>(config.mesh.n_per_direction_y)
            : static_cast<std::size_t>(config.mesh.n_per_direction);
        sites = s2d::make_sites(
            boundary,
            s2d::CartesianGridCount2D{nx_cells, ny_cells},
            s2d::SiteValidationOptions{});
    } else if (config.mesh.type == "hexagonal_square" || config.mesh.type == "voronoi_hexagonal") {
        sites = s2d::make_sites(
            boundary,
            s2d::HexagonalGrid2D{n_dir},
            s2d::SiteValidationOptions{});
    } else if (config.mesh.type == "random_square" || config.mesh.type == "voronoi_random") {
        const auto n_random = static_cast<std::size_t>(
            config.mesh.n_generators > 0 ? config.mesh.n_generators : config.mesh.n_per_direction * config.mesh.n_per_direction);
        sites = s2d::make_sites(
            boundary,
            s2d::UniformRandom2D{n_random, static_cast<std::uint32_t>(config.mesh.random_seed)},
            s2d::SiteValidationOptions{});
    } else {
        throw std::runtime_error("unknown mesh.type: " + config.mesh.type);
    }

    const auto diagram = vd2d::ClippedVoronoiBuilder2D::build(sites, boundary);
    const Index n = static_cast<Index>(diagram.volume_count());
    MeshTopology topo;
    topo.neighbor_start.reserve(static_cast<std::size_t>(n) + 1);
    topo.coord_x.assign(static_cast<std::size_t>(n), 0.0);
    topo.coord_y.assign(static_cast<std::size_t>(n), 0.0);
    topo.is_boundary_cell.assign(static_cast<std::size_t>(n), 0U);
    std::vector<Real> volumes(static_cast<std::size_t>(n), 0.0);

    for (const auto& cell : diagram.cells) {
        const auto p = static_cast<std::size_t>(cell.volume_id);
        const auto site_index = static_cast<std::size_t>(cell.site_id.value);
        const auto& site = diagram.sites[site_index];
        topo.coord_x[p] = site.point.x;
        topo.coord_y[p] = site.point.y;
        topo.is_boundary_cell[p] = static_cast<unsigned char>(cell.is_boundary_cell);
        volumes[p] = std::max<Real>(cell.area(), Real{1.0e-30});
    }

    topo.neighbor_start.push_back(0);
    for (Index p = 0; p < n; ++p) {
        const auto& cell = diagram.cells[static_cast<std::size_t>(p)];
        const auto ps = static_cast<std::size_t>(cell.volume_id);
        const Real px = topo.coord_x[ps];
        const Real py = topo.coord_y[ps];
        for (const auto& edge : cell.edges) {
            if (edge.is_boundary_edge || edge.length <= Real{0}) {
                continue;
            }
            const Index q = neighbor_volume_for_edge(diagram, cell, edge, px, py);
            if (q < 0) {
                continue;
            }
            if (q == static_cast<Index>(cell.volume_id)) {
                continue;
            }
            const auto qs = static_cast<std::size_t>(q);
            const Real dist = distance_between(
                topo.coord_x[ps], topo.coord_y[ps],
                topo.coord_x[qs], topo.coord_y[qs]);
            topo.neighbor_list.push_back(q);
            topo.interface_area.push_back(std::max<Real>(edge.length, Real{1.0e-30}));
            topo.distance.push_back(std::max(dist, Real{1.0e-30}));
        }
        topo.neighbor_start.push_back(static_cast<Index>(topo.neighbor_list.size()));
    }

    std::vector<Index> perm = reorder_mesh(config, topo, volumes);
    MeshAdapter mesh{std::move(topo), std::move(volumes)};
    mesh.set_ordering(std::move(perm));
    return mesh;
}

Real MeshAdapter::h_eff() const noexcept {
    return n_cells() > 0 ? std::sqrt(total_volume() / static_cast<Real>(n_cells())) : 0.0;
}

Real MeshAdapter::total_volume() const noexcept {
    return std::accumulate(data_.volume.begin(), data_.volume.end(), Real{0});
}

void MeshAdapter::validate() const {
    const auto n = static_cast<std::size_t>(data_.n_cells);
    if (topology_.neighbor_start.size() != n + 1) {
        throw std::runtime_error("mesh CSR neighbor_start must have n_cells + 1 entries");
    }
    if (topology_.coord_x.size() != n || topology_.coord_y.size() != n) {
        throw std::runtime_error("mesh coordinate arrays must match n_cells");
    }
    if (topology_.neighbor_list.size() != topology_.interface_area.size()
        || topology_.neighbor_list.size() != topology_.distance.size()) {
        throw std::runtime_error("mesh CSR edge arrays must have equal length");
    }
    for (const Real v : data_.volume) {
        if (v <= 0) {
            throw std::runtime_error("cell volume must be positive");
        }
    }
    for (std::size_t e = 0; e < topology_.neighbor_list.size(); ++e) {
        if (topology_.interface_area[e] <= 0 || topology_.distance[e] <= 0) {
            throw std::runtime_error("interface area and distance must be positive");
        }
    }
}

} // namespace gss
