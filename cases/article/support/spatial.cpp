// papers/PaperIII/FadeevConvergence/pap_p3_fadeev_convergence.cpp
//
// Paper III, verification against the analytical Fadeev solution.
//
// This is a PAPER program, not a study.  Its purpose is deliberately narrow:
// verify spatial convergence of the finite-volume discretisation against the
// analytical rotated Fadeev equilibrium.  Because the exact solution is known,
// the source dF/dPsi = exp(-2 Psi_F) is evaluated at the analytical field and
// kept fixed.  The resulting linear manufactured Poisson solve isolates the
// spatial operator and boundary treatment; no gradient reconstruction and no
// nonlinear iteration are involved in this verification.
//
// Outputs are written under Result/run_<UTC>/ in the program directory.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "gsmaker/analysis/error/ConvergenceOrder.hpp"
#include "gsmaker/analysis/spectral/SlepcContext.hpp"

#include "common/ConsoleTee.hpp"
#include "common/FadeevNonlinearStudyCommon.hpp"

namespace ps = study;

using Real = study::Real;
using Index = study::Index;
using DenseVector = study::DenseVector;

namespace {

struct Level {
    Index nx{0};
    Index ny{0};
};

struct Settings {
    std::vector<std::string> mesh_families{"cartesian", "hexagonal", "random"};
    std::vector<Level> levels{{48, 16}, {96, 32}, {192, 64}};
    std::uint32_t mesh_seed{12345U};
    Real alpha{static_cast<Real>(0.225)};
    Real theta_degrees{static_cast<Real>(5.7)};
    bool require_monotone_l2{true};
};

struct ErrorNorms {
    Real l2{0};
    Real linf{0};
};

struct Row {
    std::string mesh_family;
    Index nominal_nx{0};
    Index nominal_ny{0};
    Index cells{0};
    Real h{0};
    Real l2{0};
    Real linf{0};
    Index linear_iterations{0};
    Real linear_true_relative_residual{0};
    bool converged{false};
};

struct FamilySummary {
    std::string mesh_family;
    Real p_l2{std::numeric_limits<Real>::quiet_NaN()};
    bool order_available{false};
    bool monotone_l2{false};
    bool monotone_linf{false};
    bool all_converged{false};
};

[[nodiscard]] Settings load_settings(const std::string& path) {
    Settings settings;
    const YAML::Node root = YAML::LoadFile(path);

    if (const auto mesh = root["mesh"]) {
        if (mesh["families"]) {
            settings.mesh_families.clear();
            for (const auto& item : mesh["families"]) {
                settings.mesh_families.push_back(item.as<std::string>());
            }
        }
        if (mesh["random_seed"]) {
            settings.mesh_seed = mesh["random_seed"].as<std::uint32_t>();
        }
    }

    if (const auto paper = root["paper"]) {
        if (paper["alpha"]) {
            settings.alpha = static_cast<Real>(paper["alpha"].as<double>());
        }
        if (paper["theta_degrees"]) {
            settings.theta_degrees =
                static_cast<Real>(paper["theta_degrees"].as<double>());
        }
        if (paper["require_monotone_l2"]) {
            settings.require_monotone_l2 = paper["require_monotone_l2"].as<bool>();
        }
        if (paper["refinement_levels"]) {
            settings.levels.clear();
            for (const auto& node : paper["refinement_levels"]) {
                if (!node.IsSequence() || node.size() != 2U) {
                    throw std::invalid_argument(
                        "paper.refinement_levels entries must be [nx, ny]");
                }
                const auto nx = static_cast<Index>(node[0].as<int>());
                const auto ny = static_cast<Index>(node[1].as<int>());
                if (!(nx > Index{0}) || !(ny > Index{0})) {
                    throw std::invalid_argument(
                        "paper.refinement_levels must contain positive dimensions");
                }
                settings.levels.push_back(Level{nx, ny});
            }
        }
    }

    if (settings.mesh_families.empty()) {
        throw std::invalid_argument("mesh.families must not be empty");
    }
    if (settings.levels.size() != 3U) {
        throw std::invalid_argument(
            "Paper III convergence table requires exactly three refinement levels");
    }
    if (!(settings.alpha > Real{0}) ||
        !std::isfinite(static_cast<double>(settings.alpha))) {
        throw std::invalid_argument("paper.alpha must be finite and positive");
    }
    if (!std::isfinite(static_cast<double>(settings.theta_degrees))) {
        throw std::invalid_argument("paper.theta_degrees must be finite");
    }

    return settings;
}

// Analytical integrated source for the exact Fadeev equilibrium.
// The finite-volume matrix represents the negative Laplacian, so the integrated
// source entering A Psi = b_D + s is -V_i exp(-2 Psi_F,i).
[[nodiscard]] DenseVector analytical_integrated_source_at_exact(
    const ps::CaseWorkspace& workspace) {

    if (workspace.exact.size() != workspace.volumes.size()) {
        throw std::invalid_argument(
            "analytical_integrated_source_at_exact size mismatch");
    }

    DenseVector source(workspace.exact.size(), Real{0});
    for (std::size_t i = 0; i < workspace.exact.size(); ++i) {
        const double psi = static_cast<double>(workspace.exact[i]);
        const double value = std::exp(-2.0 * psi);
        if (!std::isfinite(value)) {
            throw std::runtime_error(
                "analytical Fadeev source is non-finite at the exact field");
        }
        source[i] = -workspace.volumes[i] * static_cast<Real>(value);
        if (!std::isfinite(static_cast<double>(source[i]))) {
            throw std::runtime_error(
                "integrated analytical Fadeev source is non-finite");
        }
    }
    return source;
}

struct LinearVerificationResult {
    DenseVector solution;
    Index iterations{0};
    Real true_relative_residual{0};
    bool converged{false};
    std::string reason;
};

[[nodiscard]] LinearVerificationResult solve_manufactured_poisson(
    const ps::CaseWorkspace& workspace) {

    const auto source = analytical_integrated_source_at_exact(workspace);
    if (workspace.poisson.boundary_rhs.size() != source.size()) {
        throw std::runtime_error("boundary/source size mismatch");
    }

    DenseVector rhs(source.size(), Real{0});
    for (std::size_t i = 0; i < source.size(); ++i) {
        rhs[i] = workspace.poisson.boundary_rhs[i] + source[i];
        if (!std::isfinite(static_cast<double>(rhs[i]))) {
            throw std::runtime_error("manufactured Fadeev RHS is non-finite");
        }
    }

    gss::PetscLinearSystemController linear{
        workspace.poisson.matrix, ps::direct_linear_config(true)};
    const auto solved = linear.solve(rhs);

    LinearVerificationResult result;
    result.solution = solved.solution;
    result.iterations = solved.iterations;
    result.true_relative_residual =
        static_cast<Real>(solved.true_relative_residual);
    result.converged = solved.converged();
    result.reason = solved.reason_string;
    return result;
}

[[nodiscard]] ErrorNorms solution_error(
    std::span<const Real> numerical,
    std::span<const Real> exact,
    std::span<const Real> volumes) {

    if (numerical.size() != exact.size() || numerical.size() != volumes.size()) {
        throw std::invalid_argument("solution_error size mismatch");
    }
    if (numerical.empty()) {
        throw std::invalid_argument("solution_error requires a non-empty field");
    }

    long double weighted_square_sum = 0.0L;
    long double volume_sum = 0.0L;
    long double maximum = 0.0L;

    for (std::size_t i = 0; i < numerical.size(); ++i) {
        const long double volume = static_cast<long double>(volumes[i]);
        if (!(volume > 0.0L) || !std::isfinite(static_cast<double>(volume))) {
            throw std::runtime_error("non-positive or non-finite control-volume measure");
        }

        const long double error =
            static_cast<long double>(numerical[i]) -
            static_cast<long double>(exact[i]);

        weighted_square_sum += volume * error * error;
        volume_sum += volume;
        maximum = std::max(maximum, std::abs(error));
    }

    if (!(volume_sum > 0.0L)) {
        throw std::runtime_error("invalid total control-volume measure");
    }

    return ErrorNorms{
        static_cast<Real>(std::sqrt(weighted_square_sum / volume_sum)),
        static_cast<Real>(maximum)};
}

[[nodiscard]] const char* latex_family(const std::string& family) {
    if (family == "cartesian") return "Cartesian";
    if (family == "hexagonal") return "Hexagonal";
    if (family == "random") return "Random";
    return nullptr;
}

[[nodiscard]] std::string display_family(const std::string& family) {
    if (const auto* name = latex_family(family)) return name;
    return family;
}

[[nodiscard]] std::vector<Row> rows_for_family(
    const std::vector<Row>& rows,
    const std::string& family) {

    std::vector<Row> result;
    for (const auto& row : rows) {
        if (row.mesh_family == family) result.push_back(row);
    }
    std::sort(result.begin(), result.end(), [](const Row& a, const Row& b) {
        return a.cells < b.cells;
    });
    return result;
}

[[nodiscard]] FamilySummary summarise_family(
    const std::vector<Row>& rows,
    const std::string& family) {

    const auto family_rows = rows_for_family(rows, family);
    if (family_rows.size() != 3U) {
        throw std::runtime_error(
            "expected exactly three convergence rows for mesh family " + family);
    }

    if (!(family_rows[0].cells < family_rows[1].cells &&
          family_rows[1].cells < family_rows[2].cells)) {
        throw std::runtime_error(
            "actual cell counts are not successively refined for mesh family " + family);
    }

    FamilySummary summary;
    summary.mesh_family = family;
    summary.all_converged = std::all_of(
        family_rows.begin(), family_rows.end(),
        [](const Row& row) { return row.converged; });

    if (!summary.all_converged) {
        return summary;
    }

    summary.monotone_l2 =
        family_rows[0].l2 > family_rows[1].l2 &&
        family_rows[1].l2 > family_rows[2].l2;
    summary.monotone_linf =
        family_rows[0].linf > family_rows[1].linf &&
        family_rows[1].linf > family_rows[2].linf;

    // Do not ask the nonlinear three-mesh estimator for an order when the
    // prerequisite monotone decrease is absent.  In that case the publication
    // table explicitly reports no order instead of manufacturing one.
    if (summary.monotone_l2) {
        summary.p_l2 = gss::ConvergenceOrder::estimate(
            family_rows[0].cells, family_rows[0].l2,
            family_rows[1].cells, family_rows[1].l2,
            family_rows[2].cells, family_rows[2].l2);
        summary.order_available = std::isfinite(static_cast<double>(summary.p_l2));
    }

    return summary;
}




void write_compact_latex_table(
    const std::filesystem::path& path,
    const std::vector<Row>& rows,
    const std::vector<FamilySummary>& summaries) {

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("could not create compact LaTeX table");
    }

    std::map<std::string, FamilySummary> summary_by_family;
    for (const auto& summary : summaries) {
        summary_by_family[summary.mesh_family] = summary;
    }

    out << "% Generated by pap_p3_fadeev_convergence\n"
        << "% Analytical Fadeev source evaluated at Psi_F; direct manufactured Poisson solve.\n"
        << "\\begin{table}[t]\n"
        << "\\centering\n"
        << "\\caption{Convergence of the finite-volume solution to the analytical "
           "Fadeev equilibrium. $N_{\\rm vol}$ is the actual number of generated "
           "control volumes and $p_2$ is computed from the three volume-weighted "
           "$L_2$ errors when they decrease monotonically.}\n"
        << "\\label{tab:fadeev_convergence}\n"
        << "\\begin{tabular}{lrrrr}\n"
        << "\\toprule\n"
        << "Mesh & $N_{\\rm vol}$ & $E_2$ & $E_\\infty$ & $p_2$ \\\\\n"
        << "\\midrule\n";

    out << std::scientific << std::setprecision(3);

    for (std::size_t family_index = 0; family_index < summaries.size(); ++family_index) {
        const auto& family = summaries[family_index].mesh_family;
        const auto family_rows = rows_for_family(rows, family);
        const auto& summary = summary_by_family.at(family);

        for (std::size_t i = 0; i < family_rows.size(); ++i) {
            const auto& row = family_rows[i];
            out << (i == 0U ? display_family(family) : std::string{}) << " & "
                << row.cells << " & "
                << row.l2 << " & "
                << row.linf << " & ";

            if (i == family_rows.size() - 1U && summary.order_available) {
                out << summary.p_l2;
            } else if (i == family_rows.size() - 1U) {
                out << "--";
            }
            out << " \\\\\n";
        }

        if (family_index + 1U < summaries.size()) {
            out << "\\midrule\n";
        }
    }

    out << "\\bottomrule\n"
        << "\\end{tabular}\n"
        << "\\end{table}\n";
}

void write_summary_csv(
    const std::filesystem::path& path,
    const std::vector<FamilySummary>& summaries) {

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("could not create convergence_summary.csv");
    }
    out << std::scientific << std::setprecision(17);
    out << "mesh_family,p_l2,order_available,monotone_l2,monotone_linf,"
           "all_linear_solves_converged\n";
    for (const auto& summary : summaries) {
        out << summary.mesh_family << ','
            << summary.p_l2 << ','
            << (summary.order_available ? 1 : 0) << ','
            << (summary.monotone_l2 ? 1 : 0) << ','
            << (summary.monotone_linf ? 1 : 0) << ','
            << (summary.all_converged ? 1 : 0) << '\n';
    }
}

} // namespace

int spatial_main() {
#ifndef GSMAKER_USE_SLEPC
    std::cerr << "pap_p3_fadeev_convergence requires GSMakerLib built with SLEPc\n";
    return 1;
#else
    try {
        // Existing Paper III pattern: SLEPc initialisation also initialises PETSc.
        gss::SlepcContext slepc;

        const auto settings = load_settings("config.yaml");
        const auto base_config = gss::YAMLConfig::load("config.yaml");
        const auto study_settings = ps::load_study_settings("config.yaml");

        const auto result_dir = ps::create_result_directory();
        study::ConsoleTee console(result_dir / "console.log");

        std::ofstream csv(result_dir / "fadeev_convergence.csv");
        std::ofstream manifest(result_dir / "manifest.txt");
        if (!csv || !manifest) {
            throw std::runtime_error("could not create result files");
        }

        csv << std::scientific << std::setprecision(17);
        csv << "mesh_family,nominal_nx,nominal_ny,cells,h,"
               "E2_volume_weighted,Einf,linear_iterations,linear_true_relative_residual,converged\n";

        manifest
            << "program=pap_p3_fadeev_convergence\n"
            << "paper=PaperIII\n"
            << "purpose=finite_volume_verification_against_analytical_fadeev_solution\n"
            << "source_evaluation=analytical_exp_minus_2psi_at_exact_fadeev_state\n"
            << "gradient_reconstruction=none\n"
            << "domain_aspect_ratio=3\n"
            << "mesh_size_definition=h=1/sqrt(Nvol)\n"
            << "order_metric=volume_weighted_L2\n"
            << "run_state=in_progress\n"
            << "started=" << ps::utc_run_id() << '\n';
        manifest.flush();

        const gss::Fadeev fadeev{
            settings.alpha, settings.theta_degrees, Real{1}};

        std::vector<Row> rows;
        std::size_t failures = 0U;

        std::cout
            << "\nPaper III: finite-volume verification against analytical Fadeev\n"
            << "  alpha              : " << settings.alpha << '\n'
            << "  theta [deg]        : " << settings.theta_degrees << '\n'
            << "  source             : analytical exp(-2 Psi_F), frozen at exact state\n"
            << "  nonlinear solve    : none (manufactured Poisson verification)\n"
            << "  gradient method    : none\n"
            << "  mesh families      : " << settings.mesh_families.size() << '\n'
            << "  refinement levels  : " << settings.levels.size() << "\n\n";

        for (const auto& family : settings.mesh_families) {
            std::cout << "  " << display_family(family) << '\n';

            for (const auto& level : settings.levels) {
                const auto mesh_config = ps::make_mesh_config(
                    base_config.mesh,
                    family,
                    settings.mesh_seed,
                    level.nx,
                    level.ny);

                auto workspace = ps::make_workspace(
                    mesh_config,
                    fadeev,
                    study_settings.jacobian_stencil_radius);

                const auto result = solve_manufactured_poisson(workspace);

                const auto cells = static_cast<Index>(workspace.exact.size());
                const Real h = gss::ConvergenceOrder::mesh_size(cells);

                Row row;
                row.mesh_family = family;
                row.nominal_nx = level.nx;
                row.nominal_ny = level.ny;
                row.cells = cells;
                row.h = h;
                row.linear_iterations = result.iterations;
                row.linear_true_relative_residual = result.true_relative_residual;
                row.converged = result.converged;

                if (result.converged) {
                    const auto error = solution_error(
                        result.solution,
                        workspace.exact,
                        workspace.volumes);
                    row.l2 = error.l2;
                    row.linf = error.linf;
                } else {
                    row.l2 = std::numeric_limits<Real>::quiet_NaN();
                    row.linf = std::numeric_limits<Real>::quiet_NaN();
                    ++failures;
                }

                rows.push_back(row);

                csv << family << ','
                    << level.nx << ',' << level.ny << ','
                    << row.cells << ',' << row.h << ','
                    << row.l2 << ',' << row.linf << ','
                    << row.linear_iterations << ','
                    << row.linear_true_relative_residual << ','
                    << (row.converged ? 1 : 0) << '\n';
                csv.flush();

                std::cout
                    << "    nominal=" << level.nx << 'x' << level.ny
                    << " actual_cells=" << row.cells
                    << " E2=" << std::scientific << std::setprecision(8) << row.l2
                    << " Einf=" << row.linf
                    << " KSP_it=" << row.linear_iterations
                    << " lin_rel=" << row.linear_true_relative_residual
                    << (row.converged ? " PASS" : " FAIL")
                    << '\n';
            }
            std::cout << '\n';
        }

        std::vector<FamilySummary> summaries;
        summaries.reserve(settings.mesh_families.size());

        for (const auto& family : settings.mesh_families) {
            try {
                auto summary = summarise_family(rows, family);

                if (!summary.all_converged) {
                    ++failures;
                }
                if (settings.require_monotone_l2 && !summary.monotone_l2) {
                    ++failures;
                }
                if (summary.monotone_l2 && !summary.order_available) {
                    ++failures;
                }

                std::cout
                    << "  summary " << std::left << std::setw(10) << family
                    << " p_L2=";
                if (summary.order_available) {
                    std::cout << std::scientific << std::setprecision(8)
                              << summary.p_l2;
                } else {
                    std::cout << "n/a";
                }
                std::cout
                    << " monotone_L2=" << (summary.monotone_l2 ? "yes" : "no")
                    << " monotone_Linf=" << (summary.monotone_linf ? "yes" : "no")
                    << " all_converged=" << (summary.all_converged ? "yes" : "no")
                    << '\n';

                summaries.push_back(std::move(summary));
            } catch (const std::exception& error) {
                ++failures;
                std::cerr
                    << "  convergence summary failed for " << family
                    << ": " << error.what() << '\n';

                FamilySummary summary;
                summary.mesh_family = family;
                summaries.push_back(summary);
            }
        }

        write_summary_csv(result_dir / "convergence_summary.csv", summaries);
        write_compact_latex_table(
            result_dir / "table_fadeev_convergence.tex", rows, summaries);

        manifest
            << "run_state=completed\n"
            << "validation_state=" << (failures == 0U ? "passed" : "failed") << '\n'
            << "failures=" << failures << '\n'
            << "alpha=" << settings.alpha << '\n'
            << "theta_degrees=" << settings.theta_degrees << '\n'
            << "mesh_seed=" << settings.mesh_seed << '\n'
            << "require_monotone_l2=" << (settings.require_monotone_l2 ? 1 : 0) << '\n'
            << "outputs=fadeev_convergence.csv;convergence_summary.csv;"
               "table_fadeev_convergence.tex;console.log;manifest.txt\n";

        std::cout
            << "\nResults: " << result_dir.string() << '\n'
            << "Failures: " << failures << '\n';

        return failures == 0U ? 0 : 1;

    } catch (const std::exception& error) {
        std::cerr
            << "pap_p3_fadeev_convergence failed: "
            << error.what() << '\n';
        return 1;
    }
#endif
}
