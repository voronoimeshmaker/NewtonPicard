// SPDX-License-Identifier: GPL-3.0-or-later
/** @file test_frozen_poisson.cpp
 * @brief Regression tests for boundary geometry and alpha-beta-gamma assembly.
 * @ingroup tests
 */
#include <cmath>
#include <vector>

#include <gtest/gtest.h>
#include <VoronoiMeshMaker/Boundary2D/Builders/Boundary2DBuilder.hpp>
#include <VoronoiMeshMaker/Boundary2D/Shapes/Rectangle.hpp>
#include <VoronoiMeshMaker/Sites2D/Factory/SiteFactory.hpp>
#include <VoronoiMeshMaker/Voronoi2D/Diagram/ClippedVoronoiBuilder2D.hpp>

#include "gss/assembly/FrozenPoisson2D.hpp"
#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/linear/CholeskyBackend.hpp"
#include "gss/linear/PETScLinearSystem.hpp"

namespace {

[[nodiscard]] auto make_diagram(std::size_t nx, std::size_t ny) {
    const auto domain = vmm::b2d::make_boundary(
        vmm::b2d::Rectangle(vmm::s2d::Point2{-1, -1}, 2, 2));
    const auto sites = vmm::s2d::make_sites(
        domain, vmm::s2d::CartesianGridCount2D{nx, ny},
        vmm::s2d::SiteValidationOptions{});
    vmm::vd2d::ClippedVoronoiBuildOptions2D options{};
    options.allow_parallel_cell_build = false;
    return vmm::vd2d::ClippedVoronoiBuilder2D::build(sites, domain, options);
}

} // namespace

TEST(BoundaryCondition, CoefficientsUseReadOnlyApplicationPoint) {
    const gss::BoundaryApplicationPoint2D point{7, 2, 3, 0, 1};
    const gss::DirichletBoundary2D dirichlet{
        [](const auto& p) { return p.x() + p.y(); }};
    const gss::NeumannBoundary2D neumann{
        [](const auto& p) { return p.outward_normal_y() * 4; }};
    const gss::RobinBoundary2D robin{
        2, 5, [](const auto& p) { return static_cast<gss::Real>(p.volume_id()); }};

    const auto d = dirichlet.coefficients_at(point);
    const auto n = neumann.coefficients_at(point);
    const auto r = robin.coefficients_at(point);
    EXPECT_EQ(d.alpha, 1);
    EXPECT_EQ(d.beta, 0);
    EXPECT_EQ(d.gamma, 5);
    EXPECT_EQ(n.alpha, 0);
    EXPECT_EQ(n.beta, 1);
    EXPECT_EQ(n.gamma, 4);
    EXPECT_EQ(r.alpha, 2);
    EXPECT_EQ(r.beta, 5);
    EXPECT_EQ(r.gamma, 7);
}

TEST(BoundaryGeometry, ApplicationPointIsNormalToFace) {
    const auto diagram = make_diagram(12, 10);
    const auto boundary = gss::frozen_boundary_faces(diagram);
    ASSERT_FALSE(boundary.empty());
    for (const auto& face : boundary) {
        const auto& point = face.application_point();
        const auto& cell = diagram.cells[static_cast<std::size_t>(face.volume_id())];
        const auto generator = diagram.sites[
            static_cast<std::size_t>(cell.site_id.value)].point;
        const gss::Real dx = point.x() - generator.x;
        const gss::Real dy = point.y() - generator.y;
        EXPECT_NEAR(std::hypot(dx, dy), face.perpendicular_distance(), 1e-13);
        EXPECT_NEAR(dx * point.outward_normal_y()
                        - dy * point.outward_normal_x(),
                    0, 1e-13);
    }
}

TEST(FrozenPoisson, LinearPatchAndGenericBoundaryForms) {
    ASSERT_EQ(PetscInitializeNoArguments(), 0);
    const auto diagram = make_diagram(12, 10);
    const auto mesh = gss::frozen_mesh_adapter(diagram);
    const auto geometry = gss::frozen_boundary_faces(diagram);
    const auto dirichlet = gss::evaluate_boundary_condition(
        geometry,
        gss::DirichletBoundary2D{[](const auto& point) {
            return 3 + 2 * point.x() - 5 * point.y();
        }});
    Mat matrix = gss::assemble_frozen_poisson_matrix(mesh, dirichlet);
    std::vector<gss::Real> exact(diagram.volume_count());
    std::vector<gss::Real> source(exact.size(), 0);
    for (const auto& cell : diagram.cells) {
        const auto p = diagram.sites[
            static_cast<std::size_t>(cell.site_id.value)].point;
        exact[cell.volume_id] = 3 + 2 * p.x - 5 * p.y;
    }
    const auto rhs = gss::assemble_frozen_poisson_rhs(mesh, dirichlet, source);
    Vec x = gss::PETScLinearSystem::make_vec(exact);
    Vec b = gss::PETScLinearSystem::make_vec(rhs);
    Vec residual = gss::PETScLinearSystem::make_vec(mesh.n_cells());
    ASSERT_EQ(MatMult(matrix, x, residual), 0);
    ASSERT_EQ(VecAXPY(residual, -1, b), 0);
    PetscReal norm = 0;
    ASSERT_EQ(VecNorm(residual, NORM_INFINITY, &norm), 0);
    EXPECT_LT(norm, 1e-11);

    const auto incompatible_neumann = gss::evaluate_boundary_condition(
        geometry, gss::NeumannBoundary2D{[](const auto&) { return 2; }});
    EXPECT_THROW(static_cast<void>(gss::assemble_frozen_poisson_rhs(
        mesh, incompatible_neumann, source)), std::runtime_error);
    const auto neumann = gss::evaluate_boundary_condition(
        geometry, gss::NeumannBoundary2D{[](const auto& point) {
            return 2 * point.outward_normal_y();
        }});
    const auto neumann_rhs = gss::assemble_frozen_poisson_rhs(
        mesh, neumann, source);
    std::vector<gss::Real> expected_neumann(source.size(), 0);
    for (const auto& face : geometry) {
        expected_neumann[static_cast<std::size_t>(face.volume_id())]
            += 2 * face.application_point().outward_normal_y() * face.measure();
    }
    EXPECT_EQ(neumann_rhs, expected_neumann);

    const auto robin = gss::evaluate_boundary_condition(
        geometry, gss::RobinBoundary2D{2, 5, [](const auto&) { return 7; }});
    const auto robin_rhs = gss::assemble_frozen_poisson_rhs(mesh, robin, source);
    std::vector<gss::Real> expected_robin(source.size(), 0);
    for (const auto& face : geometry) {
        expected_robin[static_cast<std::size_t>(face.volume_id())]
            += 7 * face.measure() / (2 * face.perpendicular_distance() + 5);
    }
    EXPECT_EQ(robin_rhs, expected_robin);

    const auto& first = geometry.front();
    const std::vector<gss::BoundaryFaceContribution2D> singular{
        {first,
         gss::BoundaryCoefficients2D{
             1, -first.perpendicular_distance(), 0},
         1}};
    EXPECT_THROW(
        static_cast<void>(gss::assemble_frozen_poisson_matrix(mesh, singular)),
        std::runtime_error);

    std::fill(source.begin(), source.end(), 4);
    const auto positive_laplacian =
        gss::assemble_frozen_poisson_rhs(mesh, dirichlet, source);
    for (std::size_t i = 0; i < source.size(); ++i) {
        EXPECT_NEAR(positive_laplacian[i] - rhs[i],
                    -4 * mesh.data().volume[i], 1e-13);
    }
    PetscBool symmetric = PETSC_FALSE;
    ASSERT_EQ(MatIsSymmetric(matrix, 1e-12, &symmetric), 0);
    EXPECT_TRUE(symmetric);

    {
        gss::Config solver_config;
        solver_config.solver.linear = "cholesky_sparse";
        solver_config.petsc.mat_solver_type = "mumps";
        solver_config.petsc.mumps_reorder = "pord";
        gss::CholeskyBackend solver(solver_config);
        solver.factorize(matrix);
        PC pc{};
        Mat factor{};
        PetscInt requested_ordering{};
        PetscInt effective_ordering{};
        ASSERT_EQ(KSPGetPC(solver.native_ksp(), &pc), 0);
        ASSERT_EQ(PCFactorGetMatrix(pc, &factor), 0);
        ASSERT_EQ(MatMumpsGetIcntl(factor, 7, &requested_ordering), 0);
        ASSERT_EQ(MatMumpsGetInfog(factor, 7, &effective_ordering), 0);
        EXPECT_EQ(requested_ordering, 4);
        EXPECT_EQ(effective_ordering, 4);
    }
    VecDestroy(&residual);
    VecDestroy(&b);
    VecDestroy(&x);
    MatDestroy(&matrix);
    EXPECT_EQ(PetscFinalize(), 0);
}
