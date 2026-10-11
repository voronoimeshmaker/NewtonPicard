// SPDX-License-Identifier: GPL-3.0-or-later
/** @file FrozenPoisson2D.cpp
 * @brief Shared source assembly and alpha-beta-gamma boundary operator.
 * @ingroup assembly
 */
#include "gss/assembly/FrozenPoisson2D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/linear/PETScLinearSystem.hpp"

namespace gss {

MeshAdapter frozen_mesh_adapter(
    const vmm::vd2d::ClippedVoronoiDiagram2D& diagram) {
    MeshTopology topology;
    std::vector<Real> volumes;
    topology.neighbor_start.push_back(0);
    for (const auto& cell : diagram.cells) {
        if (cell.volume_id != volumes.size()) {
            throw std::runtime_error("Non-contiguous volume ordering");
        }
        const auto p = diagram.sites[static_cast<std::size_t>(cell.site_id.value)].point;
        topology.coord_x.push_back(p.x);
        topology.coord_y.push_back(p.y);
        topology.is_boundary_cell.push_back(
            static_cast<unsigned char>(cell.is_boundary_cell));
        volumes.push_back(cell.area());
        for (const auto& edge : cell.edges) {
            if (edge.is_boundary_edge) {
                continue;
            }
            if (edge.length <= 0
                || edge.neighbour_site_id == vmm::s2d::kInvalidSiteId) {
                throw std::runtime_error("Invalid FV face metadata");
            }
            const auto& neighbour = diagram.cell(edge.neighbour_site_id);
            const auto q = diagram.sites[
                static_cast<std::size_t>(neighbour.site_id.value)].point;
            topology.neighbor_list.push_back(
                static_cast<Index>(neighbour.volume_id));
            topology.interface_area.push_back(edge.length);
            topology.distance.push_back(std::hypot(q.x - p.x, q.y - p.y));
        }
        topology.neighbor_start.push_back(
            static_cast<Index>(topology.neighbor_list.size()));
    }
    return MeshAdapter(std::move(topology), std::move(volumes));
}

std::vector<BoundaryFaceGeometry2D> frozen_boundary_faces(
    const vmm::vd2d::ClippedVoronoiDiagram2D& diagram) {
    std::vector<BoundaryFaceGeometry2D> result;
    const Real epsilon = Real{64} * std::numeric_limits<Real>::epsilon();
    for (const auto& cell : diagram.cells) {
        const Real local_scale = std::sqrt(cell.area());
        for (const auto& edge : cell.edges) {
            if (!edge.is_boundary_edge) {
                continue;
            }
            const auto& geometry = edge.boundary_condition_geometry();
            if (!geometry.valid() || !geometry.normal_matches_face()) {
                std::ostringstream message;
                message << "Invalid boundary-condition application geometry at volume "
                        << cell.volume_id << ": face_length=" << edge.length
                        << ", d_perp=" << geometry.perpendicular_distance()
                        << ", valid=" << geometry.valid()
                        << ", normal_matches_face="
                        << geometry.normal_matches_face();
                throw std::runtime_error(message.str());
            }
            if (!(edge.length > epsilon * local_scale)) {
                throw std::runtime_error("Degenerate boundary face measure");
            }
            if (!(geometry.perpendicular_distance() > epsilon * local_scale)) {
                std::ostringstream message;
                message << "Degenerate boundary perpendicular distance at volume "
                        << cell.volume_id;
                throw std::runtime_error(message.str());
            }
            const auto point = geometry.point();
            const auto normal = geometry.outward_normal();
            const auto volume_id = static_cast<Index>(cell.volume_id);
            result.emplace_back(
                volume_id,
                edge.length,
                geometry.perpendicular_distance(),
                BoundaryApplicationPoint2D{
                    volume_id, point.x, point.y, normal.x, normal.y},
                geometry.point_inside_local_face());
        }
    }
    return result;
}

BoundaryGeometryDiagnostics2D diagnose_frozen_boundary_faces(
    const vmm::vd2d::ClippedVoronoiDiagram2D& diagram) {
    const auto faces = frozen_boundary_faces(diagram);
    BoundaryGeometryDiagnostics2D result;
    if (faces.empty()) {
        return result;
    }
    result.min_measure_over_distance = std::numeric_limits<Real>::infinity();
    result.min_measure_over_sqrt_cell_area = std::numeric_limits<Real>::infinity();
    Real offset_sum = 0;
    for (const auto& cell : diagram.cells) {
        for (const auto& edge : cell.edges) {
            if (!edge.is_boundary_edge) {
                continue;
            }
            const auto& face = faces[result.total_faces++];
            result.application_points_outside_face += !face.application_point_inside_face();
            const auto& point = face.application_point();
            const auto midpoint = edge.midpoint;
            const Real offset = std::hypot(point.x() - midpoint.x, point.y() - midpoint.y)
                / std::sqrt(face.measure());
            offset_sum += offset;
            result.max_application_offset_over_sqrt_measure =
                std::max(result.max_application_offset_over_sqrt_measure, offset);
            const Real ratio = face.measure() / face.perpendicular_distance();
            result.min_measure_over_distance = std::min(result.min_measure_over_distance, ratio);
            result.max_measure_over_distance = std::max(result.max_measure_over_distance, ratio);
            result.min_measure_over_sqrt_cell_area = std::min(
                result.min_measure_over_sqrt_cell_area, face.measure() / std::sqrt(cell.area()));
        }
    }
    result.mean_application_offset_over_sqrt_measure = offset_sum / static_cast<Real>(result.total_faces);
    return result;
}

namespace {

[[nodiscard]] Real uniform_boundary_diffusivity(
    std::span<const BoundaryFaceContribution2D> boundary) {
    const Real diffusivity = boundary.empty() ? Real{1} : boundary.front().diffusivity();
    if (!(diffusivity > Real{0}) || !std::isfinite(diffusivity)) {
        throw std::runtime_error("Boundary diffusivity must be finite and positive");
    }
    for (const auto& face : boundary) {
        if (face.diffusivity() != diffusivity) {
            throw std::runtime_error("Only uniform diffusivity is supported: boundary faces have different values");
        }
    }
    return diffusivity;
}

[[nodiscard]] bool is_pure_neumann(std::span<const BoundaryFaceContribution2D> boundary) {
    return std::all_of(boundary.begin(), boundary.end(), [](const auto& face) {
        return face.coefficients().alpha == Real{0};
    });
}

[[nodiscard]] Real boundary_denominator(
    const BoundaryFaceContribution2D& face) {
    const auto& coefficients = face.coefficients();
    const Real denominator = coefficients.alpha
        * face.geometry().perpendicular_distance() + coefficients.beta;
    const Real scale = std::abs(coefficients.alpha)
        * face.geometry().perpendicular_distance()
        + std::abs(coefficients.beta);
    if (!(std::abs(denominator)
          > Real{64} * std::numeric_limits<Real>::epsilon()
                * std::max(scale, Real{1}))) {
        throw std::runtime_error("Singular boundary equation alpha*d + beta");
    }
    return denominator;
}

} // namespace

Mat assemble_frozen_poisson_matrix(
    const MeshAdapter& mesh,
    std::span<const BoundaryFaceContribution2D> boundary) {
    const Real diffusivity = uniform_boundary_diffusivity(boundary);
    if (is_pure_neumann(boundary)) {
        throw std::runtime_error("Pure Neumann problem: the solution is defined only up to a constant; anchor psi in one volume or fix its mean before assembling a nonsingular matrix");
    }
    std::vector<Real> denominators;
    denominators.reserve(boundary.size());
    for (const auto& face : boundary) {
        denominators.push_back(boundary_denominator(face));
    }
    bool has_positive_boundary_diagonal = false;
    bool has_negative_boundary_diagonal = false;
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const auto& face = boundary[i];
        const Real diagonal = face.diffusivity() * face.geometry().measure()
            * face.coefficients().alpha / denominators[i];
        has_positive_boundary_diagonal =
            has_positive_boundary_diagonal || diagonal > Real{0};
        has_negative_boundary_diagonal =
            has_negative_boundary_diagonal || diagonal < Real{0};
    }
    Mat matrix = MatrixAssembler::assemble(mesh, PETSC_COMM_SELF, diffusivity);
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        const auto& face = boundary[i];
        const auto& geometry = face.geometry();
        const Real diagonal = face.diffusivity() * geometry.measure()
            * face.coefficients().alpha / denominators[i];
        petsc_check(
            MatSetValue(matrix, geometry.volume_id(), geometry.volume_id(),
                        diagonal, ADD_VALUES),
            "add boundary diagonal");
    }
    petsc_check(MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY),
                "assemble boundary diagonal");
    petsc_check(MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY),
                "finish boundary diagonal");
    if (has_positive_boundary_diagonal && !has_negative_boundary_diagonal) {
        petsc_check(MatSetOption(matrix, MAT_SPD, PETSC_TRUE),
                    "declare anchored boundary operator SPD");
    }
    return matrix;
}

std::vector<Real> assemble_frozen_poisson_rhs(
    const MeshAdapter& mesh,
    std::span<const BoundaryFaceContribution2D> boundary,
    std::span<const Real> source) {
    if (source.size() != static_cast<std::size_t>(mesh.n_cells())) {
        throw std::runtime_error("Frozen RHS array size mismatch");
    }
    const Real diffusivity = uniform_boundary_diffusivity(boundary);
    if (is_pure_neumann(boundary)) {
        long double source_integral = 0;
        long double boundary_integral = 0;
        long double scale = 0;
        for (std::size_t i = 0; i < source.size(); ++i) {
            const long double term = static_cast<long double>(mesh.data().volume[i]) * source[i];
            source_integral += term;
            scale += std::abs(term);
        }
        for (const auto& face : boundary) {
            const long double term = static_cast<long double>(diffusivity)
                * face.geometry().measure() * face.coefficients().gamma / boundary_denominator(face);
            boundary_integral += term;
            scale += std::abs(term);
        }
        // Outward flux equals the volume source for div(k grad psi) = source.
        const long double imbalance = boundary_integral - source_integral;
        const long double tolerance = 1.0e-12L * scale;
        if (!std::isfinite(imbalance) || std::abs(imbalance) > tolerance) {
            std::ostringstream message;
            message.precision(17);
            message << "Incompatible pure Neumann data: outward_flux - volume_source="
                    << imbalance << ", tolerance=" << tolerance
                    << "; anchor psi or fix its mean only after correcting this imbalance";
            throw std::runtime_error(message.str());
        }
    }
    std::vector<Real> rhs(source.size());
    // MatrixAssembler represents minus the integrated Laplacian.
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        rhs[i] = -mesh.data().volume[i] * source[i];
    }
    for (const auto& face : boundary) {
        const auto& geometry = face.geometry();
        const Real contribution = face.diffusivity() * geometry.measure()
            * face.coefficients().gamma / boundary_denominator(face);
        rhs[static_cast<std::size_t>(geometry.volume_id())] += contribution;
    }
    return rhs;
}

} // namespace gss
