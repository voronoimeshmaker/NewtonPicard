#include "../cases/fadeev/common/FadeevDirichletRHS.hpp"
#include <gtest/gtest.h>
#include <fstream>
#include <iomanip>
#include "gss/solvers/NewtonSolver.hpp"
#include "gss/solvers/PicardSolver.hpp"
#include "gss/linear/CholeskyBackend.hpp"

namespace {
class PetscEnvironment : public ::testing::Environment {
    void SetUp() override { gss::petsc_check(PetscInitializeNoArguments(), "initialize Fadeev test"); }
    void TearDown() override { PetscFinalize(); }
};
[[maybe_unused]] auto* env = ::testing::AddGlobalTestEnvironment(new PetscEnvironment);

using fadeev_case::FadeevDirichletRHS;
}

TEST(FadeevSolvers, NewtonThenPicardRefinement) {
    gss::FadeevSolution exact;
    double previous=0;
    for(int n : {8,16,32}) {
        gss::Config c;
        c.mesh.domain={-1,1,-1,1}; c.mesh.n_generators=n*n; c.mesh.n_per_direction=n;
        c.solver.diagonal_shift=0; c.solver.max_iter=200;
        c.solver.tolerance_absolute=1.e-13; c.solver.tolerance_relative=1.e-12;
        auto mesh=gss::MeshAdapter::uniform_from_config(c);
        std::vector<double> initial(mesh.n_cells(),0);
        gss::NewtonSolver<gss::CholeskyBackend,gss::ADOLCJacobian,FadeevDirichletRHS> newton(c);
        const auto nr=newton.solve(mesh,initial);
        ASSERT_TRUE(nr.converged) << "Newton n=" << n;
        gss::PicardSolver<gss::CholeskyBackend,FadeevDirichletRHS> picard(c);
        const auto pr=picard.solve(mesh,initial);
        ASSERT_TRUE(pr.converged) << "Picard n=" << n;
        double difference=0;
        for(std::size_t i=0;i<initial.size();++i) difference=std::max(difference,std::abs(nr.psi[i]-pr.psi[i]));
        EXPECT_LT(difference,1.e-9);
        double new_error=0;
        for(int solver=0;solver<2;++solver) {
            const auto& r=solver==0?nr:pr;
            std::vector<double> error(r.psi.size());
            for(std::size_t i=0;i<error.size();++i) {
                const double x=mesh.topology().coord_x[i], y=mesh.topology().coord_y[i];
                const double analytical=exact.psi_exact(x,y);
                EXPECT_NEAR(exact.laplacian(x,y),(1.1/0.9)*std::exp(-2*analytical),1.e-13);
                error[i]=r.psi[i]-analytical;
            }
            const auto norms=gss::compute_norms(error,mesh.data().volume);
            FadeevDirichletRHS rhs(c);
            const auto b=rhs.assemble(mesh,r.psi);
            const auto a=gss::MatrixAssembler::apply_stencil(mesh,r.psi);
            double residual=0;
            for(std::size_t i=0;i<a.size();++i) residual=std::max(residual,std::abs(a[i]-b[i]));
            EXPECT_LT(residual,1.e-10);
            const double order=previous>0?std::log2(previous/norms.L2):0;
            if(previous>0) { EXPECT_GT(order,1.8); }
            if(solver==0) new_error=norms.L2;
            std::cout << n << " " << (solver==0?"Newton":"Picard") << " iterations=" << r.iterations << " L2=" << norms.L2 << " Linf=" << norms.Linf << " order=" << order << '\n';
        }
        previous=new_error;
    }
}