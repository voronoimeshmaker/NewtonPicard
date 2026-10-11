#pragma once
#include "FadeevDirichletRHS.hpp"
#include "gss/io/YAMLConfig.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace fadeev_case {
inline int run(const std::string& method, int argc, char** argv) {
    PetscInitialize(&argc,&argv,nullptr,nullptr);
    int code=0;
    try {
        int ranks; MPI_Comm_size(PETSC_COMM_WORLD,&ranks);
        if(ranks!=1) throw std::runtime_error("Run these report cases with one MPI process");
        const std::filesystem::path case_dir=FADEEV_CASE_DIR;
        const std::filesystem::path config=argc>1?std::filesystem::path(argv[1]):case_dir/"config.yaml";
        const auto output=std::filesystem::weakly_canonical(case_dir/"results");
        if(argc>2 && std::filesystem::weakly_canonical(argv[2])!=output)
            throw std::runtime_error("Outputs must remain in this case results directory");
        if(output.parent_path()!=std::filesystem::canonical(case_dir))
            throw std::runtime_error("Results directory resolves outside the case directory");
        auto c=gss::YAMLConfig::load(config);
        if(c.solver.nonlinear!=method || c.solver.diagonal_shift!=0
           || c.mesh.domain!=std::array<double,4>{-1,1,-1,1}
           || c.benchmark.fadeev_epsilon!=0.1 || c.benchmark.fadeev_psi0!=1 || c.benchmark.fadeev_k!=1
           || c.mesh.type!="uniform_internal" || c.solver.initial_guess!="zero"
           || !c.convergence_study.enabled || c.convergence_study.n_generators_sequence.empty())
            throw std::runtime_error("Configuration incompatible with this fixed Fadeev benchmark; see README");
        std::filesystem::create_directories(output);
        gss::YAMLConfig::write(c,output/"config_used.yaml");
        std::ofstream summary(output/"convergence.csv");
        summary.exceptions(std::ios::badbit|std::ios::failbit);
        summary << "n,cells,iterations,L1,L2,Linf,residual_Linf,order_L2,seconds\n" << std::setprecision(16);
        double previous=0,previous_h=0;
        gss::FadeevSolution exact;
        for(auto cells:c.convergence_study.n_generators_sequence) {
            const int n=static_cast<int>(std::sqrt(cells));
            if(n<2 || n*n!=cells || (previous_h>0 && 2.0/n>=previous_h))
                throw std::runtime_error("Use increasing square cell counts, e.g. 64,256,1024");
            c.mesh.n_generators=cells; c.mesh.n_per_direction=n; c.mesh.n_per_direction_y=n;
            auto mesh=gss::MeshAdapter::uniform_from_config(c);
            std::vector<double> initial(cells,0);
            gss::SolveResult result;
            if(method=="newton") {
                gss::NewtonSolver<gss::CholeskyBackend,gss::ADOLCJacobian,FadeevDirichletRHS> solver(c);
                result=solver.solve(mesh,initial);
            } else {
                gss::PicardSolver<gss::CholeskyBackend,FadeevDirichletRHS> solver(c);
                result=solver.solve(mesh,initial);
            }
            if(!result.converged) throw std::runtime_error("Solver did not converge at n="+std::to_string(n));
            FadeevDirichletRHS rhs(c);
            const auto b=rhs.assemble(mesh,result.psi);
            const auto a=gss::MatrixAssembler::apply_stencil(mesh,result.psi);
            std::vector<double> error(cells);
            double residual=0;
            std::ofstream field(output/("solution_"+std::to_string(n)+".csv"));
            field.exceptions(std::ios::badbit|std::ios::failbit);
            field << "x,y,psi_numerical,psi_analytical,error,residual\n" << std::setprecision(16);
            for(int i=0;i<cells;++i) {
                const auto x=mesh.topology().coord_x[i],y=mesh.topology().coord_y[i];
                const auto truth=exact.psi_exact(x,y);
                error[i]=result.psi[i]-truth;
                residual=std::max(residual,std::abs(a[i]-b[i]));
                field << x << ',' << y << ',' << result.psi[i] << ',' << truth << ',' << error[i] << ',' << a[i]-b[i] << '\n';
            }
            const auto norm=gss::compute_norms(error,mesh.data().volume);
            const double h=2.0/n, order=previous>0?std::log(previous/norm.L2)/std::log(previous_h/h):0;
            summary << n << ',' << cells << ',' << result.iterations << ',' << norm.L1 << ',' << norm.L2 << ',' << norm.Linf << ',' << residual << ',' << order << ',' << result.timing.total << '\n';
            std::ofstream history(output/("residual_history_"+std::to_string(n)+".csv"));
            history.exceptions(std::ios::badbit|std::ios::failbit);
            history << "evaluation,residual_norm\n" << std::setprecision(16);
            for(std::size_t k=0;k<result.residual_history.size();++k) history << k+1 << ',' << result.residual_history[k] << '\n';
            std::cout << method << " n=" << n << " iterations=" << result.iterations << " L2=" << norm.L2 << " order=" << order << '\n';
            if(!std::isfinite(norm.L2) || !std::isfinite(residual) || residual>1.e-10 || (previous>0 && order<1.8))
                throw std::runtime_error("Residual or refinement acceptance criterion failed");
            previous=norm.L2;previous_h=h;
        }
    } catch(const std::exception& e) { std::cerr << e.what() << '\n';code=1; }
    PetscFinalize();return code;
}
}
