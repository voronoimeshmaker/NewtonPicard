// SPDX-License-Identifier: GPL-3.0-or-later
/** @file FrozenPoisson2D.hpp
 * @brief Boundary-aware assembly for the frozen equation Laplacian(psi)=source.
 * @ingroup assembly
 */
#pragma once

#include <span>
#include <vector>

#include <petscmat.h>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiDiagram2D.hpp>

#include "gss/boundary/BoundaryCondition2D.hpp"
#include "gss/core/MeshAdapter.hpp"

namespace gss {

class BoundaryFaceGeometry2D {
public:
    BoundaryFaceGeometry2D(Index volume_id,
                           Real measure,
                           Real perpendicular_distance,
                           BoundaryApplicationPoint2D application_point,
                           bool application_point_inside_face) noexcept
        : volume_id_(volume_id), measure_(measure),
          perpendicular_distance_(perpendicular_distance),
          application_point_(application_point),
          application_point_inside_face_(application_point_inside_face) {}

    [[nodiscard]] Index volume_id() const noexcept { return volume_id_; }
    [[nodiscard]] Real measure() const noexcept { return measure_; }
    [[nodiscard]] Real perpendicular_distance() const noexcept {
        return perpendicular_distance_;
    }
    [[nodiscard]] const BoundaryApplicationPoint2D& application_point() const noexcept {
        return application_point_;
    }
    [[nodiscard]] bool application_point_inside_face() const noexcept {
        return application_point_inside_face_;
    }
    [[nodiscard]] BoundaryFaceGeometry2D with_volume_id(Index volume_id) const noexcept {
        return {
            volume_id,
            measure_,
            perpendicular_distance_,
            BoundaryApplicationPoint2D{
                volume_id,
                application_point_.x(),
                application_point_.y(),
                application_point_.outward_normal_x(),
                application_point_.outward_normal_y()},
            application_point_inside_face_};
    }

private:
    Index volume_id_{};
    Real measure_{};
    Real perpendicular_distance_{};
    BoundaryApplicationPoint2D application_point_{0, 0, 0, 0, 0};
    bool application_point_inside_face_{false};
};

class BoundaryFaceContribution2D {
public:
    BoundaryFaceContribution2D(BoundaryFaceGeometry2D geometry,
                               BoundaryCoefficients2D coefficients,
                               Real diffusivity) noexcept
        : geometry_(geometry), coefficients_(coefficients),
          diffusivity_(diffusivity) {}

    [[nodiscard]] const BoundaryFaceGeometry2D& geometry() const noexcept {
        return geometry_;
    }
    [[nodiscard]] const BoundaryCoefficients2D& coefficients() const noexcept {
        return coefficients_;
    }
    [[nodiscard]] Real diffusivity() const noexcept { return diffusivity_; }

private:
    BoundaryFaceGeometry2D geometry_;
    BoundaryCoefficients2D coefficients_;
    Real diffusivity_{};
};

struct BoundaryGeometryDiagnostics2D {
    std::size_t total_faces{};
    std::size_t application_points_outside_face{};
    Real max_application_offset_over_sqrt_measure{};
    Real mean_application_offset_over_sqrt_measure{};
    Real min_measure_over_distance{};
    Real max_measure_over_distance{};
    Real min_measure_over_sqrt_cell_area{};
};

[[nodiscard]] BoundaryGeometryDiagnostics2D diagnose_frozen_boundary_faces(
    const vmm::vd2d::ClippedVoronoiDiagram2D& diagram);

[[nodiscard]] MeshAdapter frozen_mesh_adapter(
    const vmm::vd2d::ClippedVoronoiDiagram2D& diagram);
[[nodiscard]] std::vector<BoundaryFaceGeometry2D> frozen_boundary_faces(
    const vmm::vd2d::ClippedVoronoiDiagram2D& diagram);

template<BoundaryCondition2D Condition>
[[nodiscard]] std::vector<BoundaryFaceContribution2D> evaluate_boundary_condition(
    std::span<const BoundaryFaceGeometry2D> boundary,
    const Condition& condition,
    Real diffusivity = Real{1}) {
    std::vector<BoundaryFaceContribution2D> result;
    result.reserve(boundary.size());
    for (const auto& face : boundary) {
        result.emplace_back(
            face, condition.coefficients_at(face.application_point()), diffusivity);
    }
    return result;
}

[[nodiscard]] Mat assemble_frozen_poisson_matrix(
    const MeshAdapter& mesh,
    std::span<const BoundaryFaceContribution2D> boundary);
[[nodiscard]] std::vector<Real> assemble_frozen_poisson_rhs(
    const MeshAdapter& mesh,
    std::span<const BoundaryFaceContribution2D> boundary,
    std::span<const Real> source);

} // namespace gss
