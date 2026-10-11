#pragma once
#include "gss/solvers/NewtonSolver.hpp"
#include "gss/solvers/PicardSolver.hpp"
#include "gss/linear/CholeskyBackend.hpp"

namespace fadeev_case {
// Test case: epsilon=0.1, psi0=k=1, domain [-1,1]^2.
// Laplacian(psi) = C exp(-2 psi), C=(1+epsilon)/(1-epsilon).
// The existing positive diffusion stencil represents -Laplacian.
class FadeevDirichletRHS {
    gss::FadeevSolution exact_;
public:
    explicit FadeevDirichletRHS(const gss::Config&) {}
    std::pair<double,double> boundary(const gss::MeshAdapter& mesh, gss::Index i) const {
        const double h=std::sqrt(mesh.data().volume[i]);
        const double x=mesh.topology().coord_x[i], y=mesh.topology().coord_y[i];
        double diagonal=0, load=0;
        auto face=[&](double bx,double by) { diagonal+=2; load+=2*exact_.psi_exact(bx,by); };
        if(x < -1+h*0.75) face(-1,y);
        if(x >  1-h*0.75) face(1,y);
        if(y < -1+h*0.75) face(x,-1);
        if(y >  1-h*0.75) face(x,1);
        return {diagonal,load}; // face length / half-cell distance = 2
    }
    template<typename T>
    T load(const gss::MeshAdapter& mesh, std::span<const T> psi, gss::Index i) const {
        using std::exp;
        const auto [d,b]=boundary(mesh,i);
        return b - mesh.data().volume[i]*(1.1/0.9)*exp(-2.0*psi[i]);
    }
    template<typename T>
    T value_at(const gss::MeshAdapter& mesh, std::span<const T> psi, gss::Index i) const {
        return load(mesh,psi,i) - boundary(mesh,i).first*psi[i];
    }
    std::vector<double> assemble(gss::MeshAdapter& mesh, std::span<const double> psi) const {
        std::vector<double> b(psi.size());
        for(gss::Index i=0;i<mesh.n_cells();++i) b[i]=value_at(mesh,psi,i);
        return b;
    }
    void add_matrix_terms(const gss::MeshAdapter& mesh, Mat A) const {
        PetscInt begin,end; gss::petsc_check(MatGetOwnershipRange(A,&begin,&end),"boundary ownership");
        for(PetscInt i=begin;i<end;++i)
            gss::petsc_check(MatSetValue(A,i,i,boundary(mesh,i).first,ADD_VALUES),"Dirichlet diagonal");
        gss::petsc_check(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY),"boundary assembly");
        gss::petsc_check(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY),"boundary assembly");
    }
    void assemble_into(gss::MeshAdapter& mesh, Vec psi, Vec b) const {
        const auto x=gss::PETScLinearSystem::to_std_vector(psi);
        const std::span<const double> view(x);
        PetscInt begin,end; gss::petsc_check(VecGetOwnershipRange(b,&begin,&end),"RHS ownership");
        for(PetscInt i=begin;i<end;++i)
            gss::petsc_check(VecSetValue(b,i,load(mesh,view,i),INSERT_VALUES),"Fadeev RHS");
        gss::petsc_check(VecAssemblyBegin(b),"RHS assembly");
        gss::petsc_check(VecAssemblyEnd(b),"RHS assembly");
    }
};
}
