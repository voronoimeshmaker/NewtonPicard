#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
/** @file VoronoiFaceConnectivity2D.hpp
 * @brief Transactional pairing of physical Voronoi faces.
 * @ingroup voronoi2d_diagram
 */
//==============================================================================
// Description : Reciprocal connectivity derived only from final cell faces.
// License     : GNU GPL v3
//==============================================================================

//==============================================================================
// c++ includes
//==============================================================================
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

//==============================================================================
// VoronoiMeshMaker includes
//==============================================================================
#include <VoronoiMeshMaker/Sites2D/SiteSet.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Cells/VoronoiCell2D.hpp>

VORMAKER_NAMESPACE_OPEN
VORONOI2D_NAMESPACE_OPEN

struct VoronoiFaceLengthLimits2D {
    Real absolute{0};
    // Preserve every resolvable positive face unless the caller requests a
    // stricter mesh-quality threshold. No implicit face removal is performed.
    Real relative{0};

    [[nodiscard]] Real minimum_length(Real domain_area, std::size_t cell_count) const {
        if (!std::isfinite(absolute) || absolute < Real{0} ||
            !std::isfinite(relative) || relative < Real{0} ||
            !std::isfinite(domain_area) || domain_area <= Real{0} || cell_count == 0U) {
            throw std::invalid_argument("Invalid face-length limits, domain area or cell count");
        }
        const Real result = std::max(absolute,
            relative * std::sqrt(domain_area / static_cast<Real>(cell_count)));
        if (!std::isfinite(result)) {
            throw std::invalid_argument("Non-finite minimum face length");
        }
        return result;
    }
};

/**
 * @brief Pair the final internal segments, independently of Delaunay links.
 *
 * Endpoints must coincide in reverse order within floating-point resolution.
 * Point contacts, unmatched segments and non-manifold matches are not links.
 * All faces (including boundary faces) must exceed minimum_length. Failures
 * leave every cell unchanged. This operation does not repair or remove faces.
 */
struct VoronoiFaceConnectivity2D {
    static void rebuild(std::span<VoronoiCell2D> cells,
                        const ::vmm::s2d::SiteSet& sites,
                        Real minimum_length) {
        if (!std::isfinite(minimum_length) || minimum_length < Real{0} ||
            !sites.ids_are_sequential() || cells.size() != sites.size()) {
            throw std::invalid_argument("Invalid face connectivity inputs");
        }
        std::vector<bool> seen(sites.size(), false);
        Real scale{1};
        for (const auto& cell : cells) {
            if (cell.site_id.value < 0 ||
                static_cast<std::size_t>(cell.site_id.value) >= sites.size()) {
                throw std::invalid_argument("Invalid cell site identity");
            }
            const auto id = static_cast<std::size_t>(cell.site_id.value);
            if (seen[id]) throw std::invalid_argument("Duplicate cell site identity");
            seen[id] = true;
            const auto generator = sites[id].point;
            if (!finite(generator) || cell.empty() || cell.edges.size() != cell.polygon.size() ||
                !std::isfinite(cell.signed_area()) || cell.signed_area() <= Real{0}) {
                throw std::invalid_argument("Invalid cell geometry");
            }
            for (std::size_t e = 0; e < cell.edges.size(); ++e) {
                const auto& edge = cell.edges[e];
                if (!finite(edge.a) || !finite(edge.b) ||
                    !equal(edge.a, cell.polygon[e]) ||
                    !equal(edge.b, cell.polygon[(e + 1U) % cell.polygon.size()])) {
                    throw std::invalid_argument("Face endpoints do not match the cell polygon");
                }
                const Real length = std::hypot(edge.b.x - edge.a.x, edge.b.y - edge.a.y);
                if (!std::isfinite(length) || length <= minimum_length) {
                    throw std::runtime_error("Near-null face at site " +
                        std::to_string(cell.site_id.value) + ", edge " + std::to_string(e) +
                        "; mesh rejected, no connectivity returned");
                }
                scale = std::max({scale, std::abs(edge.a.x), std::abs(edge.a.y),
                                 std::abs(edge.b.x), std::abs(edge.b.y)});
            }
        }

        const Real tolerance = Real{64} * std::numeric_limits<Real>::epsilon() * scale;
        std::vector<Reference> references;
        std::unordered_map<Bucket, std::vector<std::size_t>, BucketHash> buckets;
        std::vector<std::vector<VoronoiCellEdge2D>> staged;
        staged.reserve(cells.size());
        for (std::size_t c = 0; c < cells.size(); ++c) {
            staged.push_back(cells[c].edges);
            for (std::size_t e = 0; e < staged.back().size(); ++e) {
                auto& edge = staged.back()[e];
                edge.length = std::hypot(edge.b.x - edge.a.x, edge.b.y - edge.a.y);
                edge.midpoint = {std::midpoint(edge.a.x, edge.b.x),
                                 std::midpoint(edge.a.y, edge.b.y)};
                edge.neighbour_site_id = ::vmm::s2d::kInvalidSiteId;
                if (edge.is_boundary_edge) continue;
                if (edge.length <= Real{2} * tolerance) {
                    throw std::runtime_error("Internal face is below geometric matching resolution");
                }
                buckets[bucket(edge.midpoint, tolerance)].push_back(references.size());
                references.push_back({c, e});
            }
        }

        std::vector<std::size_t> partners(references.size(), references.size());
        for (std::size_t i = 0; i < references.size(); ++i) {
            const auto ref = references[i];
            const auto& edge = staged[ref.cell][ref.edge];
            const auto key = bucket(edge.midpoint, tolerance);
            std::size_t matches = 0;
            for (std::int64_t dx = -1; dx <= 1; ++dx) {
                for (std::int64_t dy = -1; dy <= 1; ++dy) {
                    const auto found = buckets.find({key.x + dx, key.y + dy});
                    if (found == buckets.end()) continue;
                    for (const auto j : found->second) {
                        const auto other = references[j];
                        if (other.cell == ref.cell) continue;
                        const auto& candidate = staged[other.cell][other.edge];
                        if (close(edge.a, candidate.b, tolerance) &&
                            close(edge.b, candidate.a, tolerance)) {
                            partners[i] = j;
                            ++matches;
                        }
                    }
                }
            }
            if (matches != 1U) {
                Real best = std::numeric_limits<Real>::infinity();
                std::size_t best_c=0, best_e=0;
                for(std::size_t cc=0;cc<staged.size();++cc) {
                    if(cc==ref.cell) continue;
                    for(std::size_t ee=0;ee<staged[cc].size();++ee) {
                        const auto& cand=staged[cc][ee];
                        Real d=std::max(std::hypot(edge.a.x-cand.b.x,edge.a.y-cand.b.y),std::hypot(edge.b.x-cand.a.x,edge.b.y-cand.a.y));
                        if(d<best){best=d;best_c=cc;best_e=ee;}
                    }
                }
                const auto& cand=staged[best_c][best_e];
                std::cerr<<std::setprecision(17)<<"DIAGNOSTIC tolerance="<<tolerance<<" length="<<edge.length<<" best_endpoint_distance="<<best<<" neighbour="<<cells[best_c].site_id.value<<" edge="<<best_e<<" boundary="<<cand.is_boundary_edge<<"\n";
                std::cerr<<"original "<<edge.a.x<<","<<edge.a.y<<" -> "<<edge.b.x<<","<<edge.b.y<<"\npartner "<<cand.a.x<<","<<cand.a.y<<" -> "<<cand.b.x<<","<<cand.b.y<<"\n";

                throw std::runtime_error("Physical face must have exactly one reciprocal segment: site " +
                    std::to_string(cells[ref.cell].site_id.value) + ", edge " +
                    std::to_string(ref.edge) + ", matches " + std::to_string(matches));
            }
        }
        for (std::size_t i = 0; i < references.size(); ++i) {
            if (partners[partners[i]] != i) {
                throw std::runtime_error("Non-reciprocal geometric face pairing");
            }
            const auto owner = references[i];
            const auto neighbour = references[partners[i]];
            auto& edge = staged[owner.cell][owner.edge];
            edge.neighbour_site_id = cells[neighbour.cell].site_id;
            define_internal_geometry(edge,
                sites[static_cast<std::size_t>(cells[owner.cell].site_id.value)].point,
                sites[static_cast<std::size_t>(edge.neighbour_site_id.value)].point);
        }
        for (std::size_t c = 0; c < cells.size(); ++c) cells[c].edges.swap(staged[c]);
    }

private:
    using Point2 = ::vmm::s2d::Point2;
    struct Reference { std::size_t cell; std::size_t edge; };
    struct Bucket {
        std::int64_t x;
        std::int64_t y;
        bool operator==(const Bucket&) const = default;
    };
    struct BucketHash {
        std::size_t operator()(Bucket key) const noexcept {
            const auto x = std::hash<std::int64_t>{}(key.x);
            const auto y = std::hash<std::int64_t>{}(key.y);
            return x ^ (y + std::size_t{0x9e3779b9U} + (x << 6U) + (x >> 2U));
        }
    };
    [[nodiscard]] static bool finite(Point2 p) noexcept {
        return std::isfinite(p.x) && std::isfinite(p.y);
    }
    [[nodiscard]] static bool equal(Point2 a, Point2 b) noexcept {
        return a.x == b.x && a.y == b.y;
    }
    [[nodiscard]] static bool close(Point2 a, Point2 b, Real tolerance) noexcept {
        return std::hypot(a.x - b.x, a.y - b.y) <= tolerance;
    }
    [[nodiscard]] static Bucket bucket(Point2 p, Real tolerance) noexcept {
        // The scale-derived tolerance bounds these coordinates below 2^47.
        return {static_cast<std::int64_t>(std::floor(p.x / tolerance)),
                static_cast<std::int64_t>(std::floor(p.y / tolerance))};
    }
    static void define_internal_geometry(VoronoiCellEdge2D& edge, Point2 site, Point2 neighbour) {
        const Real dx = edge.b.x - edge.a.x;
        const Real dy = edge.b.y - edge.a.y;
        const Real gx = neighbour.x - site.x;
        const Real gy = neighbour.y - site.y;
        const Real denominator = gx * dy - gy * dx;
        if (!std::isfinite(denominator) || denominator == Real{0}) {
            throw std::runtime_error("Degenerate generator/face crossing");
        }
        const Real ax = edge.a.x - site.x;
        const Real ay = edge.a.y - site.y;
        const Real along_generators = (ax * dy - ay * dx) / denominator;
        edge.crossing_parameter = (ax * gy - ay * gx) / denominator;
        edge.representative_point = {site.x + along_generators * gx,
                                     site.y + along_generators * gy};
        edge.generator_distance = std::hypot(gx, gy);
        edge.face_distance = std::hypot(edge.representative_point.x - site.x,
                                        edge.representative_point.y - site.y);
        edge.representative_distance = edge.face_distance;
        edge.representative_inside_local_edge =
            edge.crossing_parameter >= Real{0} && edge.crossing_parameter <= Real{1};
        edge.representative_valid = edge.representative_inside_local_edge &&
            along_generators >= Real{0} && along_generators <= Real{1};
        if (!finite(edge.representative_point) || !std::isfinite(edge.crossing_parameter) ||
            !std::isfinite(edge.generator_distance) || !std::isfinite(edge.face_distance)) {
            throw std::runtime_error("Non-finite internal face geometry");
        }
    }
};

VORONOI2D_NAMESPACE_CLOSE
VORMAKER_NAMESPACE_CLOSE
