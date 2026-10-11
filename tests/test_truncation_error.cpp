#include <gtest/gtest.h>

#include "gss/analysis/TruncationErrorAnalyzer.hpp"
#include "gss/assembly/RHSAssembler.hpp"

TEST(TruncationError, SolutionErrorForExactFieldIsZero) {
    gss::Config config{};
    config.mesh.n_generators = 64;
    const auto mesh = gss::MeshAdapter::uniform_from_config(config);
    gss::FadeevSolution fadeev{};
    std::vector<gss::Real> exact(static_cast<std::size_t>(mesh.n_cells()));
    for (gss::Index i = 0; i < mesh.n_cells(); ++i) {
        const auto s = static_cast<std::size_t>(i);
        exact[s] = fadeev.psi_exact(mesh.topology().coord_x[s], mesh.topology().coord_y[s]);
    }
    using Analyzer = gss::TruncationErrorAnalyzer<gss::MeshAdapter, gss::RHSAssembler<>>;
    const auto error = Analyzer::compute_solution_error(mesh, exact, fadeev);
    EXPECT_NEAR(error.e.Linf, 0.0, 1.0e-14);
}
