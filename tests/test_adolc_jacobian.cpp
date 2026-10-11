#include <gtest/gtest.h>
#include <cmath>
#include "gss/solvers/NewtonSolver.hpp"
#include "gss/solvers/SolverFactory.hpp"

namespace {
class PetscEnvironment : public ::testing::Environment {
public:
    void SetUp() override { gss::petsc_check(PetscInitializeNoArguments(), "initialize tests"); }
    void TearDown() override { PetscFinalize(); }
};
[[maybe_unused]] auto* environment = ::testing::AddGlobalTestEnvironment(new PetscEnvironment);

struct CoupledSource {
    template<typename T>
    T compute_source(gss::Index p, gss::Real, gss::Real, std::span<const T> x) const {
        using std::sin;
        const auto q = (static_cast<std::size_t>(p) + 1) % x.size();
        return x[p] * x[p] * x[p] + 0.25 * sin(x[q]);
    }
};
struct CubicSource {
    template<typename T>
    T compute_source(gss::Index p, gss::Real, gss::Real, std::span<const T> x) const {
        return 1.0 - x[p] * x[p] * x[p];
    }
};
gss::MeshAdapter mesh_for_test() {
    gss::Config c;
    c.mesh.n_generators = 16;
    c.mesh.n_per_direction = 4;
    return gss::MeshAdapter::uniform_from_config(c);
}
}

TEST(ADOLCJacobian, PrescribedSourceMatchesDiffusionMatrix) {
    auto mesh = mesh_for_test();
    gss::Config c;
    gss::RHSAssembler<> rhs(c);
    std::vector<gss::Real> x(mesh.n_cells(), 0.3);
    Mat J = gss::ADOLCJacobian::assemble(mesh, rhs, x);
    Mat A = gss::MatrixAssembler::assemble(mesh);
    gss::petsc_check(MatAXPY(J, -1.0, A, DIFFERENT_NONZERO_PATTERN), "compare Jacobian");
    PetscReal norm;
    MatNorm(J, NORM_FROBENIUS, &norm);
    EXPECT_LT(norm, 1.e-12);
    MatDestroy(&J); MatDestroy(&A);
}

TEST(ADOLCJacobian, NonlinearCoupledSourceMatchesAnalyticDerivativeAfterRetaping) {
    auto mesh = mesh_for_test();
    gss::Config c;
    gss::RHSAssembler<CoupledSource> rhs(c);
    std::vector<gss::Real> x(mesh.n_cells());
    Mat A = gss::MatrixAssembler::assemble(mesh);
    std::vector<short> before; cachedTraceTags(before);
    for (int state = 0; state < 2; ++state) {
        for (std::size_t j = 0; j < x.size(); ++j) x[j] = 0.15 + 0.07*j + 0.3*state;
        Mat J = gss::ADOLCJacobian::assemble(mesh, rhs, x);
        PetscInt start, end; MatGetOwnershipRange(J, &start, &end);
        for (PetscInt i = start; i < end; ++i) {
            for (PetscInt j = 0; j < mesh.n_cells(); ++j) {
                PetscScalar a, actual;
                MatGetValues(A, 1, &i, 1, &j, &a);
                MatGetValues(J, 1, &i, 1, &j, &actual);
                double expected = a;
                if (j == i) expected -= mesh.data().volume[i] * 3*x[i]*x[i];
                if (j == (i+1)%mesh.n_cells()) expected -= mesh.data().volume[i]*0.25*std::cos(x[j]);
                EXPECT_NEAR(actual, expected, 1.e-12) << i << ", " << j;
            }
        }
        MatDestroy(&J);
    }
    std::vector<short> after; cachedTraceTags(after);
    EXPECT_EQ(before, after);
    MatDestroy(&A);
}

TEST(ADOLCJacobian, NewtonConvergesForNonlinearCubicSource) {
    int ranks; MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
    if (ranks > 1) GTEST_SKIP() << "Tests the serial Newton solve";
    auto mesh = mesh_for_test();
    gss::Config c;
    c.solver.diagonal_shift = 0.0; // This manufactured problem is unshifted.
    c.solver.max_iter = 30;
    c.solver.tolerance_absolute = 1.e-11;
    c.solver.tolerance_relative = 1.e-10;
    gss::NewtonSolver<gss::CholeskyBackend, gss::ADOLCJacobian,
                      gss::RHSAssembler<CubicSource>> solver(c);
    std::vector<gss::Real> initial(mesh.n_cells(), 0.5);
    const auto result = solver.solve(mesh, initial);
    ASSERT_TRUE(result.converged);
    for (const auto x : result.psi) EXPECT_NEAR(x, 1.0, 1.e-8);
}

TEST(ADOLCJacobian, FactoryRejectsFiniteDifferences) {
    gss::Config c;
    c.solver.nonlinear = "newton";
    EXPECT_NO_THROW((void)gss::SolverFactory::create(c));
    c.newton.jacobian = "finite_difference";
    EXPECT_THROW((void)gss::SolverFactory::create(c), std::invalid_argument);
}
// Independent oracle: fixed three-cell geometry, without MatrixAssembler::assemble.
TEST(ADOLCJacobian, IndependentThreeCellResidualAndJacobian) {
    gss::MeshTopology t;
    t.neighbor_start = {0, 1, 3, 4};
    t.neighbor_list = {1, 0, 2, 1};
    t.interface_area = {2, 2, 3, 3};
    t.distance = {1, 1, 1, 1};
    t.coord_x = {0, 1, 2}; t.coord_y = {0, 0, 0};
    t.is_boundary_cell = {1, 0, 1};
    const std::vector<double> volumes{0.2, 0.3, 0.5};
    gss::MeshAdapter mesh(std::move(t), volumes);
    gss::Config c;
    gss::RHSAssembler<CoupledSource> rhs(c);
    const double diffusion[3][3] = {{2,-2,0},{-2,5,-3},{0,-3,3}};
    double max_error = 0;
    for (const auto& x : std::vector<std::vector<double>>{{-0.4,0,1.2},{0.8,-0.7,0.2}}) {
        const auto b = rhs.assemble(mesh,x);
        const auto applied = gss::MatrixAssembler::apply_stencil(mesh,x);
        Mat J = gss::ADOLCJacobian::assemble(mesh,rhs,x);
        PetscInt begin,end; MatGetOwnershipRange(J,&begin,&end);
        for (PetscInt i=begin;i<end;++i) {
            const auto q=(i+1)%3;
            double expected_residual=-volumes[i]*(x[i]*x[i]*x[i]+0.25*std::sin(x[q]));
            for(PetscInt j=0;j<3;++j) {
                expected_residual+=diffusion[i][j]*x[j];
                double expected=diffusion[i][j];
                if(i==j) expected-=volumes[i]*3*x[i]*x[i];
                if(j==q) expected-=volumes[i]*0.25*std::cos(x[q]);
                double actual=0;
                ASSERT_EQ(MatGetValues(J,1,&i,1,&j,&actual),0);
                max_error=std::max(max_error,std::abs(actual-expected));
                EXPECT_NEAR(actual,expected,1.e-12);
            }
            EXPECT_NEAR(applied[i]-b[i],expected_residual,1.e-12);
        }
        MatDestroy(&J);
    }
    std::cout << "Independent analytic Jacobian max error: " << max_error << '\n';
}

TEST(ADOLCJacobian, ShiftedNonlinearJacobianAddsOnlyDiagonal) {
    auto mesh = mesh_for_test();
    gss::Config c;
    gss::RHSAssembler<CoupledSource> rhs(c);
    std::vector<double> x(mesh.n_cells(), 0.7);
    Mat base = gss::ADOLCJacobian::assemble(mesh, rhs, x, 0.0);
    Mat shifted = gss::ADOLCJacobian::assemble(mesh, rhs, x, 1.75);
    gss::petsc_check(MatAXPY(shifted, -1.0, base, DIFFERENT_NONZERO_PATTERN), "subtract Jacobians");
    gss::petsc_check(MatShift(shifted, -1.75), "subtract expected shift");
    PetscReal error;
    ASSERT_EQ(MatNorm(shifted, NORM_FROBENIUS, &error), 0);
    EXPECT_LT(error, 1.e-12);
    MatDestroy(&base); MatDestroy(&shifted);
}

TEST(ADOLCJacobian, ShiftedNewtonMatchesPicardAndSatisfiesResidual) {
    int ranks; MPI_Comm_size(PETSC_COMM_WORLD, &ranks);
    if (ranks > 1) GTEST_SKIP() << "Serial solver regression";
    auto mesh = mesh_for_test();
    gss::Config c;
    c.solver.diagonal_shift = 1.75;
    c.solver.max_iter = 10;
    std::vector<double> x(mesh.n_cells(), 0.4);
    gss::NewtonSolver<gss::CholeskyBackend> newton(c);
    gss::PicardSolver<gss::CholeskyBackend> picard(c);
    const auto result = newton.solve(mesh, x);
    const auto reference = picard.solve(mesh, x);
    ASSERT_TRUE(result.converged);
    ASSERT_TRUE(reference.converged);
    EXPECT_LE(result.iterations, 2);
    gss::RHSAssembler<> rhs(c);
    const auto b = rhs.assemble(mesh, result.psi);
    const auto a = gss::MatrixAssembler::apply_stencil(mesh, result.psi);
    for (std::size_t i=0; i<x.size(); ++i) {
        EXPECT_NEAR(result.psi[i], reference.psi[i], 1.e-11);
        EXPECT_NEAR(a[i] + 1.75*result.psi[i] - b[i], 0.0, 1.e-11);
    }
}
