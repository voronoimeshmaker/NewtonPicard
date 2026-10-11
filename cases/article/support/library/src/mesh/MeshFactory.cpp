#include "gsmaker/mesh/MeshFactory.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <VoronoiMeshMaker/Boundary2D/Builders/Boundary2DBuilder.hpp>
#include <VoronoiMeshMaker/Boundary2D/Shapes/Rectangle.hpp>
#include <VoronoiMeshMaker/Sites2D/Factory/SiteFactory.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiBuilder2D.hpp>

#include "gsmaker/mesh/BoundaryFaceGeometry.hpp"

namespace gss {
namespace {

[[nodiscard]]
Real distance_between(
    Real ax,
    Real ay,
    Real bx,
    Real by) noexcept {

    return std::hypot(
        bx - ax,
        by - ay);
}



void require_natural_ordering(
    const MeshConfig& config) {

    if (!config.reordering.empty() &&
        config.reordering != "natural") {

        throw std::runtime_error(
            "Mesh reordering is no longer performed by GSMakerLib. "
            "Apply the requested ordering through VoronoiMeshMaker before constructing the solver mesh.");
    }
}


void require_valid_domain(
    const MeshConfig& config) {

    if (config.domain[1] <= config.domain[0] ||
        config.domain[3] <= config.domain[2]) {

        throw std::runtime_error(
            "mesh.domain must have xmax > xmin and ymax > ymin");
    }
}


template<typename Diagram, typename Cell, typename Edge>
[[nodiscard]]
Index neighbor_volume_for_edge(
    const Diagram& diagram,
    const Cell& cell,
    const Edge& edge,
    Real px,
    Real py) {

    constexpr Real tolerance_scale =
        Real{1.0e-8};

    const Real centre_distance_squared =
        (edge.midpoint.x - px) *
            (edge.midpoint.x - px) +
        (edge.midpoint.y - py) *
            (edge.midpoint.y - py);

    Real best_residual =
        std::numeric_limits<Real>::max();

    Index best_volume =
        Index{-1};

    for (const auto neighbour_site_id :
         cell.neighbor_ids) {

        const auto neighbour_site_index =
            static_cast<std::size_t>(
                neighbour_site_id.value);

        const auto& neighbour_site =
            diagram.sites[
                neighbour_site_index];

        const Real neighbour_distance_squared =
            (edge.midpoint.x -
             neighbour_site.point.x) *
                (edge.midpoint.x -
                 neighbour_site.point.x) +
            (edge.midpoint.y -
             neighbour_site.point.y) *
                (edge.midpoint.y -
                 neighbour_site.point.y);

        const Real residual =
            std::abs(
                centre_distance_squared -
                neighbour_distance_squared);

        if (residual < best_residual) {
            best_residual = residual;

            best_volume =
                static_cast<Index>(
                    diagram.cell(
                        neighbour_site_id)
                        .volume_id);
        }
    }

    const Real tolerance =
        tolerance_scale *
        std::max(
            Real{1},
            centre_distance_squared);

    return best_residual <= tolerance
        ? best_volume
        : Index{-1};
}

} // namespace


MeshAdapter MeshFactory::create_uniform_internal(
    const MeshConfig& config) {

    const Index n =
        config.n_generators;

    if (n <= 0) {
        throw std::runtime_error(
            "mesh.n_generators must be positive");
    }

    require_valid_domain(
        config);

    const Index nx =
        static_cast<Index>(
            std::sqrt(
                static_cast<Real>(n)));

    if (nx <= 0) {
        throw std::runtime_error(
            "invalid Cartesian mesh resolution");
    }

    const Index ny =
        (nx * nx == n)
            ? nx
            : static_cast<Index>(
                  std::ceil(
                      static_cast<Real>(n) /
                      static_cast<Real>(nx)));

    const Real xmin =
        config.domain[0];

    const Real xmax =
        config.domain[1];

    const Real ymin =
        config.domain[2];

    const Real ymax =
        config.domain[3];

    const Real dx =
        (xmax - xmin) /
        static_cast<Real>(nx);

    const Real dy =
        (ymax - ymin) /
        static_cast<Real>(ny);

    MeshTopology topology;

    topology.neighbor_start.reserve(
        static_cast<std::size_t>(n) + 1U);

    topology.coord_x.resize(
        static_cast<std::size_t>(n));

    topology.coord_y.resize(
        static_cast<std::size_t>(n));

    topology.is_boundary_cell.assign(
        static_cast<std::size_t>(n),
        0U);

    std::vector<Real> volumes(
        static_cast<std::size_t>(n),
        dx * dy);

    std::vector<BoundaryFaceGeometry>
        boundary_faces;

    boundary_faces.reserve(
        static_cast<std::size_t>(
            2 * nx +
            2 * ny));

    topology.neighbor_start.push_back(
        Index{0});

    for (Index j = 0;
         j < ny;
         ++j) {

        for (Index i = 0;
             i < nx;
             ++i) {

            const Index p =
                j * nx + i;

            if (p >= n) {
                break;
            }

            const auto ps =
                static_cast<std::size_t>(p);

            const Real cx =
                xmin +
                (static_cast<Real>(i) +
                 Real{0.5}) *
                    dx;

            const Real cy =
                ymin +
                (static_cast<Real>(j) +
                 Real{0.5}) *
                    dy;

            topology.coord_x[ps] = cx;
            topology.coord_y[ps] = cy;

            const bool left =
                i == 0;

            const bool right =
                i == nx - 1 ||
                p + 1 >= n;

            const bool bottom =
                j == 0;

            const bool top =
                j == ny - 1 ||
                p + nx >= n;

            topology.is_boundary_cell[ps] =
                static_cast<unsigned char>(
                    left ||
                    right ||
                    bottom ||
                    top);

            const auto add_neighbor =
                [&](bool present,
                    Index neighbour,
                    Real area,
                    Real distance) {

                    if (!present) {
                        return;
                    }

                    topology.neighbor_list.push_back(
                        neighbour);

                    topology.interface_area.push_back(
                        area);

                    topology.distance.push_back(
                        distance);
                };

            add_neighbor(
                i > 0,
                p - 1,
                dy,
                dx);

            add_neighbor(
                i + 1 < nx &&
                    p + 1 < n,
                p + 1,
                dy,
                dx);

            add_neighbor(
                j > 0,
                p - nx,
                dx,
                dy);

            add_neighbor(
                j + 1 < ny &&
                    p + nx < n,
                p + nx,
                dx,
                dy);

            if (left) {
                boundary_faces.push_back({
                    p,
                    xmin,
                    cy,
                    Real{-1},
                    Real{0},
                    dy,
                    dx / Real{2},
                    true});
            }

            if (right) {
                boundary_faces.push_back({
                    p,
                    xmax,
                    cy,
                    Real{1},
                    Real{0},
                    dy,
                    dx / Real{2},
                    true});
            }

            if (bottom) {
                boundary_faces.push_back({
                    p,
                    cx,
                    ymin,
                    Real{0},
                    Real{-1},
                    dx,
                    dy / Real{2},
                    true});
            }

            if (top) {
                boundary_faces.push_back({
                    p,
                    cx,
                    ymax,
                    Real{0},
                    Real{1},
                    dx,
                    dy / Real{2},
                    true});
            }

            topology.neighbor_start.push_back(
                static_cast<Index>(
                    topology.neighbor_list.size()));
        }
    }

    return MeshAdapter{
        std::move(topology),
        std::move(volumes),
        std::move(boundary_faces)};
}


MeshFactory::VoronoiDiagram2D
MeshFactory::create_voronoi_diagram(
    const MeshConfig& config) {

    require_natural_ordering(
        config);

    require_valid_domain(
        config);

    if (config.type ==
        "uniform_internal") {

        throw std::runtime_error(
            "uniform_internal is an internal test mesh and has no VoronoiMeshMaker diagram");
    }

    namespace b2d = vmm::b2d;
    namespace s2d = vmm::s2d;
    namespace vd2d = vmm::vd2d;

    const Real xmin =
        config.domain[0];

    const Real xmax =
        config.domain[1];

    const Real ymin =
        config.domain[2];

    const Real ymax =
        config.domain[3];

    const auto boundary =
        b2d::make_boundary(
            b2d::Rectangle(
                b2d::Point2{
                    xmin,
                    ymin},
                xmax - xmin,
                ymax - ymin));

    const auto n_per_direction =
        static_cast<std::size_t>(
            std::max<Index>(
                config.n_per_direction,
                Index{1}));

    s2d::SiteSet sites;

    if (config.type ==
            "cartesian_square" ||
        config.type ==
            "voronoi_uniform") {

        const auto nx_cells =
            static_cast<std::size_t>(
                std::max<Index>(
                    config.n_per_direction,
                    Index{1}));

        const auto ny_cells =
            config.n_per_direction_y > 0
                ? static_cast<std::size_t>(
                      config.n_per_direction_y)
                : nx_cells;

        sites =
            s2d::make_sites(
                boundary,
                s2d::CartesianGridCount2D{
                    nx_cells,
                    ny_cells},
                s2d::SiteValidationOptions{});

    } else if (
        config.type ==
            "hexagonal_square" ||
        config.type ==
            "voronoi_hexagonal") {

        sites =
            s2d::make_sites(
                boundary,
                s2d::HexagonalGrid2D{
                    n_per_direction},
                s2d::SiteValidationOptions{});

    } else if (
        config.type ==
            "random_square" ||
        config.type ==
            "voronoi_random") {

        const Index requested =
            config.n_generators > 0
                ? config.n_generators
                : config.n_per_direction *
                      config.n_per_direction;

        if (requested <= 0) {
            throw std::runtime_error(
                "random mesh requires a positive number of generators");
        }

        sites =
            s2d::make_sites(
                boundary,
                s2d::UniformRandom2D{
                    static_cast<std::size_t>(
                        requested),
                    static_cast<std::uint32_t>(
                        config.random_seed)},
                s2d::SiteValidationOptions{});

    } else {
        throw std::runtime_error(
            "unknown mesh.type: " +
            config.type);
    }

    return
        vd2d::ClippedVoronoiBuilder2D::build(
            sites,
            boundary);
}


MeshAdapter MeshFactory::create(
    const VoronoiDiagram2D& diagram) {

    const Index n =
        static_cast<Index>(
            diagram.volume_count());

    if (n <= 0) {
        throw std::runtime_error(
            "VoronoiMeshMaker produced an empty mesh");
    }

    MeshTopology topology;

    topology.neighbor_start.reserve(
        static_cast<std::size_t>(n) + 1U);

    topology.coord_x.assign(
        static_cast<std::size_t>(n),
        Real{0});

    topology.coord_y.assign(
        static_cast<std::size_t>(n),
        Real{0});

    topology.is_boundary_cell.assign(
        static_cast<std::size_t>(n),
        0U);

    std::vector<Real> volumes(
        static_cast<std::size_t>(n),
        Real{0});

    std::vector<BoundaryFaceGeometry>
        boundary_faces;

    for (const auto& cell :
         diagram.cells) {

        const auto cell_index =
            static_cast<std::size_t>(
                cell.volume_id);

        const auto site_index =
            static_cast<std::size_t>(
                cell.site_id.value);

        const auto& site =
            diagram.sites[
                site_index];

        topology.coord_x[cell_index] =
            site.point.x;

        topology.coord_y[cell_index] =
            site.point.y;

        topology.is_boundary_cell[cell_index] =
            static_cast<unsigned char>(
                cell.is_boundary_cell);

        const Real area =
            cell.area();

        if (!std::isfinite(
                static_cast<double>(
                    area)) ||
            !(area > Real{0})) {

            throw std::runtime_error(
                "VoronoiMeshMaker produced a non-positive or non-finite cell area");
        }

        volumes[cell_index] =
            area;
    }

    topology.neighbor_start.push_back(
        Index{0});

    for (Index p = 0;
         p < n;
         ++p) {

        const auto& cell =
            diagram.cells[
                static_cast<std::size_t>(
                    p)];

        const auto cell_index =
            static_cast<std::size_t>(
                cell.volume_id);

        const Real px =
            topology.coord_x[
                cell_index];

        const Real py =
            topology.coord_y[
                cell_index];

        for (const auto& edge :
             cell.edges) {

            if (!std::isfinite(
                    static_cast<double>(
                        edge.length)) ||
                !(edge.length > Real{0})) {

                throw std::runtime_error(
                    "VoronoiMeshMaker produced a non-positive or non-finite face measure");
            }

            if (edge.is_boundary_edge) {

                const auto& geometry =
                    edge.boundary_condition_geometry();

                if (!geometry.valid() ||
                    !geometry.normal_matches_face()) {

                    throw std::runtime_error(
                        "VoronoiMeshMaker produced invalid boundary-condition geometry");
                }

                const Real d =
                    geometry.perpendicular_distance();

                if (!std::isfinite(
                        static_cast<double>(
                            d)) ||
                    !(d > Real{0})) {

                    throw std::runtime_error(
                        "boundary perpendicular distance must be finite and positive");
                }

                const auto point =
                    geometry.point();

                const auto normal =
                    geometry.outward_normal();

                boundary_faces.push_back({
                    static_cast<Index>(
                        cell.volume_id),
                    point.x,
                    point.y,
                    normal.x,
                    normal.y,
                    edge.length,
                    d,
                    geometry.point_inside_local_face()});

                continue;
            }

            const Index neighbour =
                neighbor_volume_for_edge(
                    diagram,
                    cell,
                    edge,
                    px,
                    py);

            if (neighbour < 0 ||
                neighbour ==
                    static_cast<Index>(
                        cell.volume_id)) {

                throw std::runtime_error(
                    "failed to identify the physical neighbour of an internal Voronoi face");
            }

            const auto neighbour_index =
                static_cast<std::size_t>(
                    neighbour);

            const Real distance =
                distance_between(
                    topology.coord_x[
                        cell_index],
                    topology.coord_y[
                        cell_index],
                    topology.coord_x[
                        neighbour_index],
                    topology.coord_y[
                        neighbour_index]);

            if (!std::isfinite(
                    static_cast<double>(
                        distance)) ||
                !(distance > Real{0})) {

                throw std::runtime_error(
                    "internal generator distance must be finite and positive");
            }

            topology.neighbor_list.push_back(
                neighbour);

            topology.interface_area.push_back(
                edge.length);

            topology.distance.push_back(
                distance);
        }

        topology.neighbor_start.push_back(
            static_cast<Index>(
                topology.neighbor_list.size()));
    }

    return MeshAdapter{
        std::move(topology),
        std::move(volumes),
        std::move(boundary_faces)};
}


MeshAdapter MeshFactory::create(
    const MeshConfig& config) {

    require_natural_ordering(
        config);

    if (config.type ==
        "uniform_internal") {

        return create_uniform_internal(
            config);
    }

    const auto diagram =
        create_voronoi_diagram(
            config);

    return create(
        diagram);
}

} // namespace gss
