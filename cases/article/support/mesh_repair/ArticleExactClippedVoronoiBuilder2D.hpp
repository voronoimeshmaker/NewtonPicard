#pragma once
//==============================================================================
// Name        : ArticleExactClippedVoronoiBuilder2D.hpp
// Project     : VoronoiMeshMaker (VMM)
// Module      : Voronoi2D / Diagram
// Author      : Joao Flavio Vieira de Vasconcellos
// Version     : 1.4
// Description : Builds complete clipped 2D Voronoi diagrams.
// License     : GNU GPL v3
//==============================================================================

/**
 * @file ArticleExactClippedVoronoiBuilder2D.hpp
 * @brief Builds all clipped cells for a Boundary2D / SiteSet pair.
 *
 * The clipped cell polygon remains unchanged. This builder also fills
 * per-edge metadata so downstream numerical methods can distinguish the actual
 * clipped edge from the representative point on Boundary2D used by boundary
 * formulae.
 *
 * Delaunay connectivity and physical Voronoi-face connectivity are kept
 * distinct. Delaunay neighbours remain available through the cell-level
 * Delaunay neighbourhood, while each final internal Voronoi edge stores the
 * site on the opposite side through neighbour_site_id.
 *
 * @ingroup voronoi2d_diagram
 */

//==============================================================================
// c++ includes
//==============================================================================
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

//==============================================================================
// TBB includes
//==============================================================================
#include <tbb/parallel_for.h>

//==============================================================================
// VoronoiMeshMaker includes
//==============================================================================
#include <VoronoiMeshMaker/Boundary2D/Boundary2DData.hpp>
#include <VoronoiMeshMaker/Boundary2D/Queries/Boundary2DContains.hpp>
#include <VoronoiMeshMaker/Boundary2D/Queries/Boundary2DDistance.hpp>
#include <VoronoiMeshMaker/ErrorHandling/CoreErrors.h>
#include <VoronoiMeshMaker/ErrorHandling/ErrorManager.h>
#include <VoronoiMeshMaker/ErrorHandling/Macros.h>
#include <VoronoiMeshMaker/Sites2D/SiteSet.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Cells/BoundaryConditionPoint2D.hpp>
#include "ArticleExactVoronoiCellBuilder2D.hpp"
#include <VoronoiMeshMaker/Voronoi2D/Clipping/BoundaryShortEdgeCollapse2D.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Clipping/ClippingWorkspace2D.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Delaunay/DelaunayBuilder2D.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Delaunay/DelaunaySiteIndex.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiDiagram2D.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/VoronoiFaceConnectivity2D.hpp>

VORMAKER_NAMESPACE_OPEN
VORONOI2D_NAMESPACE_OPEN

struct ArticleExactClippedVoronoiBuildOptions2D {
    DelaunayBuildOptions2D    delaunay{};
    ArticleExactVoronoiCellBuildOptions2D cells{};
    bool allow_parallel_cell_build{true};
    Real min_boundary_edge_length{Real{0}};
    BoundaryShortEdgePolicy2D boundary_short_edge_policy{
        BoundaryShortEdgePolicy2D::Keep};
    VoronoiFaceLengthLimits2D face_length_limits{};
};

struct ArticleExactClippedVoronoiBuilder2D {
    [[nodiscard]] static ClippedVoronoiDiagram2D build(
        const ::vmm::s2d::SiteSet&          sites,
        const ::vmm::b2d::Boundary2DData&   boundary,
        const ArticleExactClippedVoronoiBuildOptions2D& options = {})
    {
        validate_inputs_or_throw(
            sites,
            boundary,
            options);

        // ---------------------------------------------------------------------
        // Build the Delaunay triangulation used to construct the Voronoi cells.
        //
        // Delaunay connectivity is intentionally retained as a distinct graph.
        // It must not be replaced by physical Voronoi-face connectivity.
        // ---------------------------------------------------------------------
        auto delaunay =
            DelaunayBuilder2D::build(
                sites,
                options.delaunay);

        auto index =
            DelaunaySiteIndex::from(
                delaunay);

        // ---------------------------------------------------------------------
        // Initialise the clipped Voronoi diagram.
        // ---------------------------------------------------------------------
        ClippedVoronoiDiagram2D diagram;

        diagram.boundary = boundary;
        diagram.sites    = sites;
        diagram.delaunay = delaunay;
        diagram.cells.resize(
            sites.size());

        std::vector<std::size_t> idx(
            sites.size());

        std::iota(
            idx.begin(),
            idx.end(),
            std::size_t{0});

        const auto& cell_opts =
            options.cells;

        // ---------------------------------------------------------------------
        // Build one clipped Voronoi cell.
        //
        // This phase may run in parallel. The cell builder uses Delaunay
        // neighbours as construction candidates for the Voronoi halfplanes.
        // ---------------------------------------------------------------------
        auto build_one_cell =
            [&](std::size_t i) {

                thread_local ClippingWorkspace2D ws;

                auto& cell =
                    diagram.cells[i];

                cell =
                    ArticleExactVoronoiCellBuilder2D::build_unchecked(
                        sites,
                        boundary,
                        delaunay,
                        index,
                        sites[i].id,
                        ws,
                        cell_opts);

                // -------------------------------------------------------------
                // Apply the configured short-boundary-edge policy before the
                // final edge metadata are constructed.
                // -------------------------------------------------------------
                collapse_short_boundary_edges(
                    cell.polygon,
                    boundary,
                    BoundaryShortEdgeCollapseOptions2D{
                        options.min_boundary_edge_length,
                        options.boundary_short_edge_policy});
            };

        // ---------------------------------------------------------------------
        // Build all clipped cell polygons.
        // ---------------------------------------------------------------------
        if (options.allow_parallel_cell_build) {

            tbb::parallel_for(std::size_t{0}, sites.size(), build_one_cell);

        } else {

            std::for_each(
                idx.begin(),
                idx.end(),
                build_one_cell);
        }

        // ---------------------------------------------------------------------
        // Build final edge metadata.
        //
        // Keep this phase outside the execution-policy algorithm because the
        // metadata validation below may throw exceptions. An exception inside
        // std::execution::par_unseq would terminate the program.
        // ---------------------------------------------------------------------
        for (const auto i : idx) {

            auto& cell =
                diagram.cells[i];

            populate_cell_edges(
                cell,
                sites[i].point,
                boundary);

            cell.is_boundary_cell =
                cell_has_boundary_edge(
                    cell);

            // -----------------------------------------------------------------
            // Boundary projection metadata.
            // -----------------------------------------------------------------
            if (cell.is_boundary_cell) {

                const auto projection =
                    ::vmm::b2d::project_to_boundary(
                        sites[i].point,
                        boundary);

                cell.boundary_projection_point =
                    projection.point;

                cell.boundary_projection_distance =
                    projection.distance;

                cell.boundary_projection_valid =
                    projection.valid;

                cell.boundary_projection_inside_cell =
                    point_in_or_on_polygon(
                        projection.point,
                        cell.polygon);

            } else {

                cell.boundary_projection_point = {};

                cell.boundary_projection_distance =
                    Real{0};

                cell.boundary_projection_valid =
                    false;

                cell.boundary_projection_inside_cell =
                    false;
            }
        }

        // ---------------------------------------------------------------------
        // Rebuild diagram indices after all cells and their edge metadata are
        // complete.
        //
        // IMPORTANT:
        // Do not overwrite the Delaunay neighbour graph here. Physical
        // connectivity is already represented by the final internal edges and
        // their neighbour_site_id values.
        // ---------------------------------------------------------------------
        VoronoiFaceConnectivity2D::rebuild(diagram.cells, sites,
            options.face_length_limits.minimum_length(domain_area(boundary), diagram.cells.size()));
        diagram.rebuild_indices();

        // Log only accepted meshes, on the caller thread. The geometry remains
        // unchanged and the existing logger controls retention of warnings.
        const auto config = ::vmm::error::Config::get();
        if (config && config->min_severity <= ::vmm::error::Severity::Warning) {
            for (const auto& cell : diagram.cells) {
                for (std::size_t e = 0; e < cell.edges.size(); ++e) {
                    const auto& edge = cell.edges[e];
                    if (edge.is_boundary_edge &&
                        !edge.boundary_condition_geometry().point_inside_local_face()) {
                        ::vmm::error::ErrorRecord warning;
                        warning.severity = ::vmm::error::Severity::Warning;
                        warning.message = "VMM_BOUNDARY_PROJECTION_OUTSIDE_FACE volume="
                            + std::to_string(cell.volume_id) + " site="
                            + std::to_string(cell.site_id.value) + " face="
                            + std::to_string(e)
                            + "; unbounded normal projection retained";
                        ::vmm::error::ErrorManager::log(std::move(warning));
                    }
                }
            }
        }

        return diagram;
    }

private:
    [[nodiscard]] static Real domain_area(const ::vmm::b2d::Boundary2DData& boundary) {
        long double area = 0;
        for (::vmm::b2d::Index r = 0; r < boundary.ring_count(); ++r) {
            const auto ring = boundary.ring(r);
            long double twice_area = 0;
            const auto origin = ring.front();
            for (std::size_t i = 0; i < ring.size(); ++i) {
                const auto a = ring[i];
                const auto b = ring[(i + 1U) % ring.size()];
                twice_area += (static_cast<long double>(a.x) - origin.x) *
                                  (static_cast<long double>(b.y) - origin.y) -
                              (static_cast<long double>(b.x) - origin.x) *
                                  (static_cast<long double>(a.y) - origin.y);
            }
            const auto magnitude = std::abs(twice_area) / 2;
            area += boundary.kinds[static_cast<std::size_t>(r)] == ::vmm::b2d::LoopKind::Hole
                ? -magnitude : magnitude;
        }
        return static_cast<Real>(area);
    }

    [[nodiscard]] static bool cell_has_boundary_edge(
        const VoronoiCell2D& cell) noexcept
    {
        return std::any_of(
            cell.edges.begin(),
            cell.edges.end(),
            [](const VoronoiCellEdge2D& edge) noexcept {
                return edge.is_boundary_edge;
            });
    }

    [[nodiscard]] static bool find_parent_boundary_segment(
        ::vmm::s2d::Point2 a,
        ::vmm::s2d::Point2 b,
        const ::vmm::b2d::Boundary2DData& boundary,
        ::vmm::b2d::Index& ring_index,
        std::size_t& edge_index,
        ::vmm::s2d::Point2& segment_a,
        ::vmm::s2d::Point2& segment_b) noexcept
    {
        constexpr Real eps =
            ::vmm::constants::kEpsilon;

        for (::vmm::b2d::Index r = 0;
             r < boundary.ring_count();
             ++r) {

            const auto ring =
                boundary.ring(r);

            if (ring.size() < 2U) {
                continue;
            }

            for (std::size_t e = 0;
                 e < ring.size();
                 ++e) {

                const auto& p =
                    ring[e];

                const auto& q =
                    ring[(e + 1U) % ring.size()];

                if (::vmm::b2d::point_on_segment(
                        a,
                        p,
                        q,
                        eps) &&
                    ::vmm::b2d::point_on_segment(
                        b,
                        p,
                        q,
                        eps)) {

                    ring_index =
                        r;

                    edge_index =
                        e;

                    segment_a =
                        p;

                    segment_b =
                        q;

                    return true;
                }
            }
        }

        ring_index =
            ::vmm::b2d::kInvalid;

        edge_index =
            0U;

        segment_a = {};
        segment_b = {};

        return false;
    }

    static void populate_cell_edges(
        VoronoiCell2D& cell,
        ::vmm::s2d::Point2 site,
        const ::vmm::b2d::Boundary2DData& boundary)
    {
        cell.edges.clear();

        const std::size_t n =
            cell.polygon.size();

        if (n < 2U) {
            return;
        }

        cell.edges.reserve(
            n);

        for (std::size_t i = 0;
             i < n;
             ++i) {

            const auto& a =
                cell.polygon[i];

            const auto& b =
                cell.polygon[(i + 1U) % n];

            VoronoiCellEdge2D edge;

            edge.a =
                a;

            edge.b =
                b;

            edge.midpoint =
                ::vmm::s2d::Point2{
                    (a.x + b.x) * Real{0.5},
                    (a.y + b.y) * Real{0.5}
                };

            const Real dx =
                b.x - a.x;

            const Real dy =
                b.y - a.y;

            edge.length =
                std::sqrt(
                    dx * dx +
                    dy * dy);

            ::vmm::b2d::Index ring_id =
                ::vmm::b2d::kInvalid;

            std::size_t boundary_edge_id =
                0U;

            ::vmm::s2d::Point2 boundary_a{};
            ::vmm::s2d::Point2 boundary_b{};

            // -----------------------------------------------------------------
            // Boundary face.
            // -----------------------------------------------------------------
            if (find_parent_boundary_segment(
                    a,
                    b,
                    boundary,
                    ring_id,
                    boundary_edge_id,
                    boundary_a,
                    boundary_b)) {

                const auto projection =
                    ::vmm::b2d::project_to_segment(
                        site,
                        boundary_a,
                        boundary_b,
                        ring_id,
                        boundary_edge_id);

                edge.is_boundary_edge =
                    true;

                edge.representative_point =
                    projection.point;

                edge.representative_distance =
                    projection.distance;

                edge.representative_valid =
                    projection.valid;

                edge.boundary_ring_index =
                    projection.ring_index;

                edge.boundary_edge_index =
                    projection.edge_index;

                edge.boundary_segment_parameter =
                    projection.segment_parameter;

                edge.representative_inside_local_edge =
                    ::vmm::b2d::point_on_segment(
                        projection.point,
                        a,
                        b,
                        ::vmm::constants::kEpsilon);

                edge.define_boundary_condition_geometry(
                    site,
                    boundary_a,
                    boundary_b);
                if (!edge.boundary_condition_geometry().valid()) {
                    throw std::runtime_error("Invalid boundary application geometry at site "
                        + std::to_string(cell.site_id.value) + ", face " + std::to_string(i));
                }
            }

            // Internal connectivity is assigned only after all final polygons exist.

            cell.edges.push_back(
                edge);
        }
    }

    [[nodiscard]] static bool point_in_or_on_polygon(
        ::vmm::s2d::Point2 point,
        const std::vector<::vmm::s2d::Point2>& polygon) noexcept
    {
        if (polygon.size() < 3U) {
            return false;
        }

        bool has_positive =
            false;

        bool has_negative =
            false;

        constexpr Real tol =
            static_cast<Real>(
                1.0e-10);

        for (std::size_t i = 0;
             i < polygon.size();
             ++i) {

            const auto& a =
                polygon[i];

            const auto& b =
                polygon[
                    (i + 1U) %
                    polygon.size()];

            const Real cross =
                (b.x - a.x) *
                    (point.y - a.y) -
                (b.y - a.y) *
                    (point.x - a.x);

            if (cross > tol) {

                has_positive =
                    true;

            } else if (cross < -tol) {

                has_negative =
                    true;
            }

            if (has_positive &&
                has_negative) {

                return false;
            }
        }

        return true;
    }

    static void validate_inputs_or_throw(
        const ::vmm::s2d::SiteSet&          sites,
        const ::vmm::b2d::Boundary2DData&   boundary,
        const ArticleExactClippedVoronoiBuildOptions2D& options)
    {
        if (sites.size() < 3U) {

            VMM_THROW(
                ::vmm::error::CoreErr::InvalidArgument,
                {
                    {
                        "where",
                        "ArticleExactClippedVoronoiBuilder2D"
                    },
                    {
                        "reason",
                        "at_least_three_sites_required"
                    },
                    {
                        "count",
                        std::to_string(
                            sites.size())
                    }
                });
        }

        if (!sites.ids_are_sequential()) {

            VMM_THROW(
                ::vmm::error::CoreErr::InvalidArgument,
                {
                    {
                        "where",
                        "ArticleExactClippedVoronoiBuilder2D"
                    },
                    {
                        "reason",
                        "site_ids_must_be_sequential"
                    }
                });
        }

        if (!boundary.invariant_ok()) {

            VMM_THROW(
                ::vmm::error::CoreErr::InvalidArgument,
                {
                    {
                        "where",
                        "ArticleExactClippedVoronoiBuilder2D"
                    },
                    {
                        "reason",
                        "invalid_boundary_invariant"
                    }
                });
        }

        validate_boundary_short_edge_options_or_throw(
            BoundaryShortEdgeCollapseOptions2D{
                options.min_boundary_edge_length,
                options.boundary_short_edge_policy});
        ArticleExactVoronoiCellBuilder2D::validate_domain(boundary);
        for (const auto& site : sites) {
            if (!std::isfinite(site.point.x) || !std::isfinite(site.point.y) ||
                !::vmm::b2d::contains(site.point, boundary)) {
                throw std::invalid_argument("Sites must be finite and inside the domain");
            }
        }
        static_cast<void>(options.face_length_limits.minimum_length(domain_area(boundary), sites.size()));
    }
};

VORONOI2D_NAMESPACE_CLOSE
VORMAKER_NAMESPACE_CLOSE
