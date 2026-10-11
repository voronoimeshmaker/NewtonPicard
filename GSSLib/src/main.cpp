#include <cmath>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <variant>

#include <petscsys.h>

#include "gss/assembly/MatrixAssembler.hpp"
#include "gss/benchmark/FadeevSolution.hpp"
#include "gss/core/MeshAdapter.hpp"
#include "gss/io/ResultsWriter.hpp"
#include "gss/io/YAMLConfig.hpp"
#include "gss/linear/CholeskyBackend.hpp"
#include "gss/linear/PETScLinearSystem.hpp"
#include "gss/solvers/SolverFactory.hpp"

namespace {

[[nodiscard]] gss::SolveResult solve_constant_linear_system(
    const gss::Config& config,
    gss::MeshAdapter& mesh) {
    gss::SolveResult result{};
    result.psi.assign(static_cast<std::size_t>(mesh.n_cells()), 0.0);

    Mat A = gss::MatrixAssembler::assemble(mesh);
    if (config.solver.diagonal_shift != 0.0) {
        gss::petsc_check(MatShift(A, config.solver.diagonal_shift), "shift linear smoke matrix");
    }

    std::vector<gss::Real> rhs(
        static_cast<std::size_t>(mesh.n_cells()),
        config.solver.rhs_constant);
    Vec b = gss::PETScLinearSystem::make_vec(rhs);
    Vec x = gss::PETScLinearSystem::make_vec(result.psi);

    gss::CholeskyBackend backend(config);
    backend.factorize(A);
    backend.solve(b, x);

    result.psi = gss::PETScLinearSystem::to_std_vector(x);
    result.converged = true;
    result.iterations = 1;

    gss::PETScLinearSystem::destroy(b);
    gss::PETScLinearSystem::destroy(x);
    gss::PETScLinearSystem::destroy(A);
    mesh.data().psi = result.psi;
    return result;
}

} // namespace

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int exit_code = 0;
    try {
        const std::filesystem::path config_path = argc > 1 ? argv[1] : "config/default.yaml";
        auto config = gss::YAMLConfig::load(config_path);
        auto mesh = gss::MeshAdapter::from_config(config);
        gss::FadeevSolution fadeev{
            config.benchmark.fadeev_epsilon,
            config.benchmark.fadeev_psi0,
            config.benchmark.fadeev_k};

        gss::SolveResult result{};
        if (config.solver.nonlinear == "linear_constant") {
            result = solve_constant_linear_system(config, mesh);
        } else {
            std::vector<gss::Real> initial(static_cast<std::size_t>(mesh.n_cells()), 0.0);
            if (config.solver.initial_guess == "fadeev_perturbed") {
                const auto& topo = mesh.topology();
                for (gss::Index i = 0; i < mesh.n_cells(); ++i) {
                    const auto s = static_cast<std::size_t>(i);
                    initial[s] = fadeev.psi_exact(topo.coord_x[s], topo.coord_y[s])
                        + config.solver.perturbation_amplitude * std::sin(topo.coord_x[s] + topo.coord_y[s]);
                }
            }
            auto solver = gss::SolverFactory::create(config);
            result = std::visit([&](auto& s) { return s.solve(mesh, initial); }, solver);
        }

        if (!result.converged) {
            std::cerr << "Solver did not converge within the iteration limit\n";
            exit_code = 2;
        }
        if (config.output.write_solution) {
            gss::ResultsWriter::write_csv("solution.csv", mesh, result.psi);
        }
        if (config.output.timing_report) {
            std::cout << "Mesh type: " << config.mesh.type
                      << " n_per_direction=" << config.mesh.n_per_direction
                      << " cells=" << mesh.n_cells()
                      << " h_eff=" << mesh.h_eff() << '\n';
            std::cout << "Converged: " << std::boolalpha << result.converged
                      << " iterations=" << result.iterations
                      << " total_time=" << result.timing.total << "s\n";
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        PetscFinalize();
        return 1;
    }
    PetscFinalize();
    return exit_code;
}
