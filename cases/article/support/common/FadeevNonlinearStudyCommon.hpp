#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>

#include <VoroGradToolkit/Gradients/Common/WLSCore2D.hpp>
#include <VoroGradToolkit/Gradients/GreenGauss/ModifiedGreenGauss2D.hpp>
#include <VoroGradToolkit/Gradients/Polynomial/WLSQR2D.hpp>
#include <VoroGradToolkit/Mesh/GradientWorkspace2D.hpp>

#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/assembly/PoissonAssembler.hpp"
#include "gsmaker/boundary/BoundaryApplication.hpp"
#include "gsmaker/core/Config.hpp"
#include "gsmaker/io/YAMLConfig.hpp"
#include "gsmaker/linear_system/PetscLinearConfig.hpp"
#include "gsmaker/linear_system/PetscLinearSystemController.hpp"
#include "gsmaker/mesh/MeshFactory.hpp"
#include "gsmaker/nonlinear/newton/Newton.hpp"
#include "gsmaker/nonlinear/picard/Picard.hpp"
#include "gsmaker/semilinear/GradientRatioSource.hpp"
#include "gsmaker/semilinear/StationaryPointRecovery.hpp"
#include "gsmaker/verification/analytical/Fadeev.hpp"

#include "SparseJacobian.hpp"

#ifdef GSMAKER_USE_SLEPC
#include "gsmaker/analysis/spectral/SlepcEigenConfig.hpp"
#include "gsmaker/analysis/spectral/SlepcEigenResult.hpp"
#include "gsmaker/analysis/spectral/SlepcEigenSolverController.hpp"
#endif

namespace study {

using Real = gss::Real;
using Index = gss::Index;
using DenseVector = std::vector<Real>;
using DenseMatrix = std::vector<Real>;

/// Beyond this cell count a dense n-by-n object is refused rather than
/// attempted.  Forming one at 21760 cells would need 3.8 GB and an O(n^3)
/// factorisation; the sparse path exists precisely so that no campaign has to.
inline constexpr std::size_t dense_cell_limit = 2000U;

struct StudySettings {
    Real alpha{static_cast<Real>(gss::Fadeev::default_alpha)};
    Real theta_degrees{static_cast<Real>(gss::Fadeev::default_theta_degrees)};
    std::uint32_t random_mesh_seed{12345U};
    Real finite_difference_scale{Real{1}};

    // CHANGED: the multiplier c of the Fadeev free function.  It is the only
    // knob in this whole family that moves the Picard spectral radius by a
    // useful amount: rho(c) = c rho(1) exactly, while alpha correlates with
    // rho at -0.51 and the domain size cannot move it at all.  At c = 1 the
    // manufactured term is identically zero and the campaign is unchanged, so
    // this defaults to 1 and every published number is reproduced.
    Real nonlinearity_scale{Real{1}};

    PetscInt requested_eigenpairs{6};
    PetscReal eps_tolerance{1.0e-11};
    PetscInt eps_max_iterations{2000};

    // CHANGED: the Jacobian stencil radius is now explicit.  The measured
    // pattern of B is the face-neighbour graph (jacobian_stencil_radius:
    // dmax(B) = 1, outside2 = 0, on three families, two gradient methods and
    // three resolutions), so 1 is the justified default.  Widen to 2 when a
    // reconstruction operator changes, and verify with verify_sparse_pattern().
    std::size_t jacobian_stencil_radius{1U};

    // CHANGED: repairing the source at stationary points of Psi alters the
    // source term, so it is opt-in and must be recorded in the manifest.  A
    // campaign that enables it is not comparable with one that does not; the
    // difference between the two is the measurement of the effect.
    bool repair_stationary_points{false};
    std::size_t recovery_sweeps{5U};

    // CHANGED: the configured initial guess is now reachable from run_picard
    // and run_newton.  It defaults to false so that every published campaign
    // number is reproduced unchanged; a study that wants to sweep the initial
    // guess sets it and the manifest records the choice.
    bool start_from_configured_guess{false};

    // CHANGED: the asymptotic factor was measured from the last four residual
    // ratios, which is too short when the sequence stops near the rounding
    // floor or when the leading mode is not clearly dominant.  The length is
    // now a parameter; 4 remains the default so existing numbers are
    // reproduced, and a campaign that wants a defensible tail raises it.
    std::size_t observed_tail_ratios{4U};
};

[[nodiscard]] inline StudySettings load_study_settings(const std::string& path) {
    StudySettings result;
    const YAML::Node root = YAML::LoadFile(path);
    const YAML::Node node = root["study"];
    if (!node) return result;

    if (node["alpha"]) result.alpha = static_cast<Real>(node["alpha"].as<double>());
    if (node["theta_degrees"]) result.theta_degrees = static_cast<Real>(node["theta_degrees"].as<double>());
    if (node["random_mesh_seed"]) result.random_mesh_seed = node["random_mesh_seed"].as<std::uint32_t>();
    if (node["finite_difference_scale"]) result.finite_difference_scale = static_cast<Real>(node["finite_difference_scale"].as<double>());
    if (node["nonlinearity_scale"]) result.nonlinearity_scale = static_cast<Real>(node["nonlinearity_scale"].as<double>());
    if (node["requested_eigenpairs"]) result.requested_eigenpairs = static_cast<PetscInt>(node["requested_eigenpairs"].as<int>());
    if (node["eps_tolerance"]) result.eps_tolerance = static_cast<PetscReal>(node["eps_tolerance"].as<double>());
    if (node["eps_max_iterations"]) result.eps_max_iterations = static_cast<PetscInt>(node["eps_max_iterations"].as<int>());

    // CHANGED: new keys, all defaulted so existing configurations keep working.
    if (node["jacobian_stencil_radius"]) {
        result.jacobian_stencil_radius =
            static_cast<std::size_t>(node["jacobian_stencil_radius"].as<int>());
    }
    if (node["repair_stationary_points"]) {
        result.repair_stationary_points = node["repair_stationary_points"].as<bool>();
    }
    if (node["recovery_sweeps"]) {
        result.recovery_sweeps = static_cast<std::size_t>(node["recovery_sweeps"].as<int>());
    }
    if (node["start_from_configured_guess"]) {
        result.start_from_configured_guess = node["start_from_configured_guess"].as<bool>();
    }
    if (node["observed_tail_ratios"]) {
        result.observed_tail_ratios =
            static_cast<std::size_t>(node["observed_tail_ratios"].as<int>());
    }

    if (!(result.alpha > Real{0}) || !std::isfinite(static_cast<double>(result.alpha))) {
        throw std::invalid_argument("study alpha must be finite and strictly positive");
    }
    if (!std::isfinite(static_cast<double>(result.theta_degrees))) {
        throw std::invalid_argument("study theta_degrees must be finite");
    }
    if (!(result.finite_difference_scale > Real{0}) ||
        !std::isfinite(static_cast<double>(result.finite_difference_scale))) {
        throw std::invalid_argument("study finite_difference_scale must be finite and positive");
    }
    // Only finiteness.  Zero and negative scales are admitted on purpose; see
    // the constructor of gss::Fadeev for why.
    if (!std::isfinite(static_cast<double>(result.nonlinearity_scale))) {
        throw std::invalid_argument("study nonlinearity_scale must be finite");
    }
    if (result.jacobian_stencil_radius == 0U || result.jacobian_stencil_radius > 3U) {
        throw std::invalid_argument("study jacobian_stencil_radius must lie between 1 and 3");
    }
    if (result.observed_tail_ratios < 2U) {
        throw std::invalid_argument("study observed_tail_ratios must be at least 2");
    }
    return result;
}

/// The Fadeev field a study should build from its own settings.
///
/// Constructing gss::Fadeev directly from settings.alpha and
/// settings.theta_degrees is what every study did, and it silently drops the
/// nonlinearity scale.  Going through this function is what keeps the field,
/// the manufactured source and the workspace at the same c.
[[nodiscard]] inline gss::Fadeev make_fadeev(const StudySettings& settings) {
    return gss::Fadeev{
        settings.alpha, settings.theta_degrees, settings.nonlinearity_scale};
}

/// Refuse to run a study that reads the scale but does not implement it.
///
/// Six studies build their source with gss::GradientRatioSource::compute
/// directly instead of going through evaluate_integrated_source, so they never
/// add the manufactured term q = (1 - c) exp(-2 Psi*):
///
///   hexagonal_generator_audit, picard_fadeev, residual_assembly,
///   residual_convergence, residual_localisation, solver_choice_fadeev
///
/// At the time of writing none of them can be affected: each has its own
/// settings loader, none calls load_study_settings, and none of their
/// configuration files carries the key, so their Fadeev is always built at
/// c = 1.  This function exists for the edit that changes that.  A study that
/// gains the shared loader and inherits the key would otherwise accept
/// nonlinearity_scale in its configuration and silently ignore it, which is the
/// same defect the Fase 6 audit confirmed twice under F05 and F06: a declared
/// capability that does not control execution.
///
/// Call it immediately after loading the settings, before any mesh is built.
inline void require_unit_nonlinearity_scale(
    const StudySettings& settings,
    std::string_view study_name) {
    if (settings.nonlinearity_scale != Real{1}) {
        throw std::invalid_argument(
            std::string{"study "} + std::string{study_name} +
            " builds its source without the manufactured term and is only "
            "valid at nonlinearity_scale = 1; remove the key or use a study "
            "that goes through evaluate_integrated_source");
    }
}

[[nodiscard]] inline std::string utc_run_id() {
    const auto now = std::chrono::system_clock::now();
    const std::time_t raw = std::chrono::system_clock::to_time_t(now);
    std::tm value{};
#if defined(_WIN32)
    gmtime_s(&value, &raw);
#else
    gmtime_r(&raw, &value);
#endif
    std::ostringstream out;
    out << "run_" << std::put_time(&value, "%Y%m%dT%H%M%SZ");
    return out.str();
}

// CHANGED: two runs started inside the same second used to share a directory.
// utc_run_id() has one-second resolution and create_directories() succeeds on a
// directory that already exists, so the second run wrote its CSV files over the
// first without any error.  A suffix is appended until the name is free, which
// makes the collision visible in the directory listing instead of silent.
[[nodiscard]] inline std::filesystem::path create_result_directory() {
    const std::filesystem::path root{"Result"};
    const auto stamp = utc_run_id();

    std::filesystem::path run = root / stamp;
    for (int suffix = 1; std::filesystem::exists(run) && suffix < 1000; ++suffix) {
        run = root / (stamp + "_" + std::to_string(suffix));
    }
    if (std::filesystem::exists(run)) {
        throw std::runtime_error("could not find a free run directory under Result/");
    }

    std::filesystem::create_directories(run);
    std::ofstream latest(root / "LATEST.txt");
    if (!latest) throw std::runtime_error("could not create Result/LATEST.txt");
    latest << run.string() << '\n';
    return run;
}

[[nodiscard]] inline Real norm2(std::span<const Real> values) {
    long double sum = 0.0L;
    for (const Real value : values) {
        const long double v = static_cast<long double>(value);
        sum += v * v;
    }
    return static_cast<Real>(std::sqrt(sum));
}

[[nodiscard]] inline Real relative_difference(
    std::span<const Real> a,
    std::span<const Real> b) {
    if (a.size() != b.size()) throw std::invalid_argument("relative_difference size mismatch");
    DenseVector delta(a.size(), Real{0});
    for (std::size_t i = 0; i < a.size(); ++i) delta[i] = a[i] - b[i];
    return norm2(delta) / std::max(Real{1}, norm2(b));
}

[[nodiscard]] inline std::size_t dense_index(
    std::size_t n,
    std::size_t row,
    std::size_t column) noexcept {
    return row * n + column;
}

// CHANGED: the resolution is now a parameter.  The previous version hardcoded
// n_per_direction = 8, n_per_direction_y = 6 and n_generators = 48, silently
// overriding whatever the caller asked for.
[[nodiscard]] inline gss::MeshConfig make_mesh_config(
    const gss::MeshConfig& base,
    std::string_view family,
    std::uint32_t seed,
    Index nx = 8,
    Index ny = 6) {
    auto config = base;
    config.n_per_direction = nx;
    config.n_per_direction_y = ny;
    config.n_generators = nx * ny;
    config.random_seed = static_cast<Index>(seed);
    config.hilbert_reorder = false;
    config.reordering = "natural";

    if (family == "cartesian") config.type = "cartesian_square";
    else if (family == "hexagonal") config.type = "hexagonal_square";
    else if (family == "random") config.type = "random_square";
    else throw std::invalid_argument("unknown mesh family");
    return config;
}

struct CaseWorkspace {
    gss::MeshFactory::VoronoiDiagram2D diagram;
    gss::MeshAdapter mesh;
    vgt::GradientWorkspace2D gradients;
    gss::PoissonSystem poisson;
    DenseVector volumes;
    DenseVector exact;

    // CHANGED: the sparsity pattern and its colouring depend only on the mesh,
    // so they are built once with the workspace and reused by every Jacobian
    // assembly.  This is the boundary with Article II: mesh-invariant
    // preparation belongs before the loop, not inside it.
    sparse::Pattern jacobian_pattern;

    // CHANGED: the scale the exact field above was built at.  The manufactured
    // source q = (1 - c) exp(-2 Psi*) is a function of that field, so a source
    // evaluated with a Fadeev at one c against a workspace populated at another
    // is arithmetic for an equation nobody posed.  It produces finite numbers
    // and no error, which is why it is checked rather than trusted.
    Real nonlinearity_scale{Real{1}};

    CaseWorkspace(
        gss::MeshFactory::VoronoiDiagram2D diagram_in,
        gss::MeshAdapter mesh_in,
        gss::PoissonSystem poisson_in,
        DenseVector volumes_in,
        DenseVector exact_in,
        std::size_t stencil_radius,
        Real nonlinearity_scale_in = Real{1})
        : diagram(std::move(diagram_in)),
          mesh(std::move(mesh_in)),
          gradients(diagram),
          poisson(std::move(poisson_in)),
          volumes(std::move(volumes_in)),
          exact(std::move(exact_in)),
          jacobian_pattern(
              sparse::build_pattern(mesh.topology(), stencil_radius)),
          nonlinearity_scale(nonlinearity_scale_in) {}

    /// Record a new scale after the exact field has been repopulated at it.
    ///
    /// A many-query study that rebuilds the field per draw calls this with the
    /// scale of the Fadeev it just used, in the same place it refreshes the
    /// Dirichlet data.  Calling it without repopulating the field defeats the
    /// cross-check rather than satisfying it.
    void set_nonlinearity_scale(Real value) noexcept {
        nonlinearity_scale = value;
    }
};

[[nodiscard]] inline CaseWorkspace make_workspace(
    const gss::MeshConfig& mesh_config,
    const gss::Fadeev& fadeev,
    std::size_t stencil_radius = 1U) {
    auto diagram = gss::MeshFactory::create_voronoi_diagram(mesh_config);
    auto mesh = gss::MeshFactory::create(diagram);

    const auto n = static_cast<std::size_t>(mesh.n_cells());
    DenseVector exact(n, Real{0});
    const auto& topology = mesh.topology();
    for (std::size_t i = 0; i < n; ++i) {
        exact[i] = fadeev.psi(topology.coord_x[i], topology.coord_y[i]);
    }

    const gss::BoundaryCondition condition = gss::Dirichlet{
        [&fadeev](const gss::BoundaryPoint2D& point) {
            return fadeev.psi(point.x, point.y);
        }};
    const auto boundary = gss::BoundaryApplication::evaluate(
        mesh.geometry().boundary_faces(), condition);
    auto poisson = gss::PoissonAssembler::assemble(mesh, boundary);

    const auto volume_span = mesh.geometry().cell_volume();
    DenseVector volumes(volume_span.begin(), volume_span.end());

    // gradients is constructed only after diagram is owned by CaseWorkspace;
    // this avoids the moved-diagram lifetime bug found in the first Picard gate.
    return CaseWorkspace{
        std::move(diagram), std::move(mesh), std::move(poisson),
        std::move(volumes), std::move(exact), stencil_radius,
        fadeev.nonlinearity_scale()};
}

[[nodiscard]] inline DenseVector free_function_values(
    std::span<const Real> psi,
    const gss::Fadeev& fadeev) {
    DenseVector result(psi.size(), Real{0});
    for (std::size_t i = 0; i < psi.size(); ++i) {
        result[i] = fadeev.free_function(psi[i]);
    }
    return result;
}

// CHANGED: diagnostics the manifest needs, returned alongside the source.
struct SourceDiagnostics {
    std::size_t regularised_cells{};
    std::size_t recovered_cells{};
    std::size_t unrecoverable_cells{};
    Real gradient_scale{};
    Real regularisation_threshold{};
    Real smallest_relative_gradient{};
};

template<class GradientMethod>
[[nodiscard]] DenseVector evaluate_integrated_source(
    const CaseWorkspace& workspace,
    std::span<const Real> psi,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    bool repair_stationary_points = false,
    std::size_t recovery_sweeps = 5U,
    SourceDiagnostics* diagnostics = nullptr) {

    if (psi.size() != workspace.volumes.size() ||
        psi.size() != workspace.exact.size()) {
        throw std::invalid_argument(
            "evaluate_integrated_source: the field, the cell volumes and the "
            "exact solution must all have the same length");
    }

    // CHANGED: the manufactured source below is built from workspace.exact, so
    // the field in the workspace and the Fadeev being evaluated have to be the
    // same problem.  Both scales come from one configuration value and are
    // copied rather than computed, so exact equality is the right test; a
    // tolerance here would only hide the mistake it is meant to catch.
    if (fadeev.nonlinearity_scale() != workspace.nonlinearity_scale) {
        throw std::invalid_argument(
            "evaluate_integrated_source: the Fadeev nonlinearity scale does "
            "not match the scale the workspace field was built at");
    }

    const auto F = free_function_values(psi, fadeev);
    const auto ratio = gss::GradientRatioSource::compute(
        workspace.gradients, psi, F, method, regularisation_delta);

    std::span<const Real> derivative{ratio.values};

    // CHANGED: where grad Psi vanishes the ratio is 0/0 and the regularised
    // value is not a derivative.  On the 211-cell hexagonal mesh the cell at
    // the Fadeev stationary point carries dF/dPsi = 0.64 and receives 0.
    gss::StationaryPointRecoveryResult repaired;

    if (repair_stationary_points && ratio.regularised_cells > 0U) {
        repaired = gss::StationaryPointRecovery::apply(
            workspace.mesh.topology(), ratio, psi, recovery_sweeps);
        derivative = std::span<const Real>{repaired.values};
    }

    if (diagnostics != nullptr) {
        diagnostics->regularised_cells = ratio.regularised_cells;
        diagnostics->recovered_cells = repaired.recovered_cells;
        diagnostics->unrecoverable_cells = repaired.unrecoverable_cells;
        diagnostics->gradient_scale = ratio.gradient_scale;
        diagnostics->regularisation_threshold = ratio.regularisation_threshold;
        diagnostics->smallest_relative_gradient = ratio.smallest_relative_gradient;
    }

    DenseVector source(psi.size(), Real{0});

    /*
     * CHANGED: the manufactured source q = (1 - c) exp(-2 Psi*) is integrated
     * with the same -V factor as the reconstructed derivative, because it
     * enters the same side of the same equation.  With A psi = s(psi) + b and
     * s = -V dF/dPsi, the discrete statement is Lap(Psi) = c exp(-2 Psi) + q,
     * whose solution is Psi* for every c.
     *
     * The c = 1 case is written as its own loop rather than as a branch inside
     * the general one.  At c = 1 the density is exactly zero, so the two are
     * mathematically the same, but only the separate loop is the same sequence
     * of floating-point operations as before this parameter existed.  That is
     * what makes a run at c = 1 a regression rather than an approximation of
     * one.  It also keeps the extra exp() per cell out of the campaign that
     * does not need it.
     */
    if (fadeev.nonlinearity_scale() == Real{1}) {
        for (std::size_t i = 0; i < psi.size(); ++i) {
            source[i] = -workspace.volumes[i] * derivative[i];
        }
        return source;
    }

    for (std::size_t i = 0; i < psi.size(); ++i) {
        source[i] = -workspace.volumes[i] *
            (derivative[i] +
             fadeev.manufactured_source_density(workspace.exact[i]));
    }
    return source;
}

template<class GradientMethod>
[[nodiscard]] DenseVector residual(
    const CaseWorkspace& workspace,
    std::span<const Real> psi,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    bool repair_stationary_points = false,
    std::size_t recovery_sweeps = 5U) {
    auto result = workspace.poisson.matrix.multiply(psi);
    const auto source = evaluate_integrated_source(
        workspace, psi, fadeev, method, regularisation_delta,
        repair_stationary_points, recovery_sweeps);
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] -= source[i] + workspace.poisson.boundary_rhs[i];
    }
    return result;
}

[[nodiscard]] inline Real component_central_step(
    Real psi_value,
    Real scale = Real{1}) {
    // Delegates so that the dense and coloured assemblies cannot drift apart.
    return sparse::central_step(psi_value, scale);
}

// CHANGED: the coloured sparse assembly replaces the one-column-at-a-time
// dense construction.  Cost falls from 2n source evaluations to 2 * colours,
// and the matrix from n^2 entries to a few per row.  The dense version is kept
// below as dense_B_reference() for verification on small meshes.
template<class GradientMethod>
[[nodiscard]] std::vector<Real> sparse_B_values(
    const CaseWorkspace& workspace,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    Real finite_difference_scale = Real{1},
    bool repair_stationary_points = false,
    std::size_t recovery_sweeps = 5U) {

    auto source = [&](std::span<const Real> psi) {
        return evaluate_integrated_source(
            workspace, psi, fadeev, method, regularisation_delta,
            repair_stationary_points, recovery_sweeps);
    };

    return sparse::assemble_values(
        workspace.jacobian_pattern,
        state,
        source,
        finite_difference_scale);
}

template<class GradientMethod>
[[nodiscard]] DenseMatrix dense_B_reference(
    const CaseWorkspace& workspace,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    Real finite_difference_scale = Real{1}) {

    if (state.size() > dense_cell_limit) {
        throw std::runtime_error(
            "dense_B_reference refused: the dense path is a verification tool "
            "for small meshes; use sparse_B_values for a campaign");
    }

    auto source = [&](std::span<const Real> psi) {
        return evaluate_integrated_source(
            workspace, psi, fadeev, method, regularisation_delta);
    };

    return sparse::assemble_dense_reference(state, source, finite_difference_scale);
}

/// Retained name for callers that still expect the dense object.
template<class GradientMethod>
[[nodiscard]] DenseMatrix finite_difference_B(
    const CaseWorkspace& workspace,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    Real finite_difference_scale = Real{1}) {
    return dense_B_reference(
        workspace, state, fadeev, method,
        regularisation_delta, finite_difference_scale);
}

/// Confirm on this mesh that the assumed stencil radius is wide enough and the
/// colouring attributes the differences correctly.  Run once per mesh family
/// and gradient operator; a pattern that is too narrow is silently wrong.
template<class GradientMethod>
[[nodiscard]] sparse::PatternComparison verify_sparse_pattern(
    const CaseWorkspace& workspace,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    Real finite_difference_scale = Real{1}) {

    const auto values = sparse_B_values(
        workspace, state, fadeev, method,
        regularisation_delta, finite_difference_scale);

    const auto dense = dense_B_reference(
        workspace, state, fadeev, method,
        regularisation_delta, finite_difference_scale);

    return sparse::compare_to_dense(workspace.jacobian_pattern, values, dense);
}

[[nodiscard]] inline DenseMatrix csr_to_dense(const gss::CsrMatrix& matrix) {
    const auto n = static_cast<std::size_t>(matrix.rows());
    if (n > dense_cell_limit) {
        throw std::runtime_error("csr_to_dense refused: mesh too large for a dense copy");
    }
    DenseMatrix dense(n * n, Real{0});
    const auto rows = matrix.row_start();
    const auto columns = matrix.column();
    const auto values = matrix.value();
    for (std::size_t row = 0; row < n; ++row) {
        const auto begin = static_cast<std::size_t>(rows[row]);
        const auto end = static_cast<std::size_t>(rows[row + 1U]);
        for (std::size_t k = begin; k < end; ++k) {
            dense[dense_index(n, row, static_cast<std::size_t>(columns[k]))] = values[k];
        }
    }
    return dense;
}

[[nodiscard]] inline gss::CsrMatrix dense_to_full_csr(
    const DenseMatrix& dense,
    std::size_t n) {
    if (dense.size() != n * n) throw std::invalid_argument("dense_to_full_csr size mismatch");
    if (n > dense_cell_limit) {
        throw std::runtime_error("dense_to_full_csr refused: mesh too large for a full CSR");
    }
    std::vector<Index> row_start(n + 1U, Index{0});
    std::vector<Index> columns;
    DenseVector values;
    columns.reserve(n * n);
    values.reserve(n * n);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            columns.push_back(static_cast<Index>(col));
            values.push_back(dense[dense_index(n, row, col)]);
        }
        row_start[row + 1U] = static_cast<Index>(values.size());
    }
    return gss::CsrMatrix{
        static_cast<Index>(n), std::move(row_start), std::move(columns), std::move(values)};
}

// CHANGED: J is now assembled in the union of the patterns of A and B, so it
// carries a few entries per row instead of n.
template<class GradientMethod>
[[nodiscard]] gss::CsrMatrix finite_difference_jacobian(
    const CaseWorkspace& workspace,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real regularisation_delta,
    Real finite_difference_scale = Real{1},
    bool repair_stationary_points = false,
    std::size_t recovery_sweeps = 5U) {

    const auto values = sparse_B_values(
        workspace, state, fadeev, method, regularisation_delta,
        finite_difference_scale, repair_stationary_points, recovery_sweeps);

    return sparse::subtract_from(
        workspace.poisson.matrix, workspace.jacobian_pattern, values);
}

// CHANGED: the factorisation package is now reachable.  PetscLinearConfig has
// carried mat_solver_type all along and PetscLinearSystemController::
// configure_ksp() already calls PCFactorSetMatSolverType when it is non-empty
// and the preconditioner is a factorisation, so nothing in the library needed
// to change; this function was clearing the field unconditionally, which is
// why the mumps_reorder and mumps_icntl_14 keys in ten configuration files
// have never had any effect.
//
// The default is still the empty string, which is PETSc's own sequential LU,
// so every existing caller is bit-for-bit unchanged.  A study that wants the
// threshold pivoting of MUMPS asks for it by name and records that it did.
[[nodiscard]] inline gss::PetscLinearConfig direct_linear_config(
    bool symmetric_spd = false,
    std::string_view factor_package = std::string_view{}) {
    gss::PetscLinearConfig config;
    config.ksp_type = "preonly";
    config.pc_type = "lu";
    config.mat_solver_type.clear();
    if (!factor_package.empty()) {
        config.mat_solver_type.assign(factor_package.data(), factor_package.size());
    }
    config.relative_tolerance = 1.0e-13;
    config.absolute_tolerance = 1.0e-15;
    config.max_iterations = 50;
    config.set_from_options = false;
    config.matrix_symmetric_hint = symmetric_spd;
    config.matrix_spd_hint = symmetric_spd;
    return config;
}

[[nodiscard]] inline gss::PicardConfig picard_config_from(
    const gss::Config& base,
    Real omega = Real{1}) {
    auto config = gss::PicardConfig::from_config(base);
    config.omega = omega;
    config.validate();
    return config;
}

[[nodiscard]] inline gss::NewtonSolverConfig newton_config_from(
    const gss::Config& base,
    Real damping = Real{1}) {
    gss::NewtonSolverConfig config;
    config.max_iterations = base.solver.max_iter;
    config.relative_residual_tolerance = base.solver.tolerance_relative;
    config.absolute_residual_tolerance = base.solver.tolerance_absolute;
    config.damping = damping;
    config.validate();
    return config;
}

// CHANGED: the configured initial guess is now reachable.  The run_picard and
// run_newton entry points below still start from workspace.exact, because
// changing their default would silently move every published campaign number;
// a study that wants the configured guess calls this and passes the result to
// the _from_state variants.
[[nodiscard]] inline DenseVector configured_initial_state(
    const gss::Config& base,
    const CaseWorkspace& workspace) {
    DenseVector state = workspace.exact;

    const Real amplitude = base.solver.perturbation_amplitude;

    if (!(amplitude > Real{0})) {
        return state;
    }
    if (base.solver.perturbation_type == "none") {
        return state;
    }

    for (std::size_t i = 0; i < state.size(); ++i) {
        const double phase = 0.731 * static_cast<double>(i + 1U);
        state[i] += amplitude * static_cast<Real>(std::sin(phase));
    }
    return state;
}

// CHANGED: both linearisations now build the source through one object.
//
// WHY.  run_picard used to call evaluate_integrated_source() with the repair
// flags at their defaults while run_newton passed settings.repair_stationary_
// points.  With repair enabled the two methods solved different problems and
// the comparison between them was meaningless, silently.  Routing both through
// SourceEvaluator makes that divergence impossible to reintroduce.
template<class GradientMethod>
struct SourceEvaluator {
    const CaseWorkspace* workspace{nullptr};
    const gss::Fadeev* fadeev{nullptr};
    const GradientMethod* method{nullptr};
    Real regularisation_delta{0};
    bool repair_stationary_points{false};
    std::size_t recovery_sweeps{5U};

    [[nodiscard]] DenseVector operator()(std::span<const Real> psi) const {
        return evaluate_integrated_source(
            *workspace, psi, *fadeev, *method, regularisation_delta,
            repair_stationary_points, recovery_sweeps);
    }
};

template<class GradientMethod>
[[nodiscard]] inline SourceEvaluator<GradientMethod> make_source_evaluator(
    const CaseWorkspace& workspace,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings) {
    return SourceEvaluator<GradientMethod>{
        &workspace,
        &fadeev,
        &method,
        base.solver.source_regularisation_delta,
        settings.repair_stationary_points,
        settings.recovery_sweeps};
}

/// The state a run starts from, chosen by the study rather than by the helper.
[[nodiscard]] inline DenseVector run_initial_state(
    const CaseWorkspace& workspace,
    const gss::Config& base,
    const StudySettings& settings) {
    return settings.start_from_configured_guess
        ? configured_initial_state(base, workspace)
        : workspace.exact;
}

// CHANGED: A is invariant across the queries, and its factorisation is the
// premise of the cost argument in the manuscript, yet every call below built a
// PetscLinearSystemController of its own and paid for the factorisation again.
// Iteration counts, spectra and norms are unaffected, because refactorising the
// same matrix gives the same answer; every published wall-clock figure that
// goes through run_picard is not.
//
// The controller is now a parameter.  A study that solves once per case keeps
// calling the forms below and behaves exactly as before; a many-query study
// builds one controller per mesh family, outside the loop, and passes it in.
// That is the regime the manuscript describes.
[[nodiscard]] inline gss::PetscLinearSystemController make_picard_controller(
    const CaseWorkspace& workspace) {
    return gss::PetscLinearSystemController{
        workspace.poisson.matrix, direct_linear_config(true)};
}

/// Entry point that reuses a factorisation the caller already owns.
template<class GradientMethod>
[[nodiscard]] gss::PicardResult run_picard(
    const CaseWorkspace& workspace,
    gss::PetscLinearSystemController& linear,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real omega,
    auto&& observer) {
    const auto source = make_source_evaluator(workspace, fadeev, method, base, settings);
    const auto initial = run_initial_state(workspace, base, settings);
    return gss::Picard::solve(
        workspace.poisson.matrix,
        workspace.poisson.boundary_rhs,
        initial,
        source,
        linear,
        picard_config_from(base, omega),
        std::forward<decltype(observer)>(observer));
}

// CHANGED: settings-aware entry point, with the same parameter order as
// run_newton.  A study that enables the stationary-point repair or the
// configured initial guess must call this form for both methods.
//
// It builds a controller per call, which is what every published campaign did.
// Keeping that behaviour here means no existing number moves by accident; a
// study that wants the reuse asks for it by passing a controller.
template<class GradientMethod>
[[nodiscard]] gss::PicardResult run_picard(
    const CaseWorkspace& workspace,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real omega,
    auto&& observer) {
    auto linear = make_picard_controller(workspace);
    return run_picard(
        workspace, linear, fadeev, method, base, settings, omega,
        std::forward<decltype(observer)>(observer));
}

/// Legacy entry point.  Equivalent to the settings-aware form with a default
/// StudySettings, which is exactly the behaviour every published campaign ran:
/// no stationary-point repair, started from the sampled analytic solution.
template<class GradientMethod>
[[nodiscard]] gss::PicardResult run_picard(
    const CaseWorkspace& workspace,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    Real omega,
    auto&& observer) {
    return run_picard(
        workspace, fadeev, method, base, StudySettings{}, omega,
        std::forward<decltype(observer)>(observer));
}

template<class GradientMethod>
[[nodiscard]] gss::NewtonResult run_newton(
    const CaseWorkspace& workspace,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real damping,
    auto&& observer) {
    const Real delta = base.solver.source_regularisation_delta;
    const auto source = make_source_evaluator(workspace, fadeev, method, base, settings);
    auto residual_evaluator = [&](std::span<const Real> psi) {
        auto value = workspace.poisson.matrix.multiply(psi);
        const auto integrated = source(psi);
        for (std::size_t i = 0; i < value.size(); ++i) {
            value[i] -= integrated[i] + workspace.poisson.boundary_rhs[i];
        }
        return value;
    };
    auto jacobian_evaluator = [&](std::span<const Real> psi) {
        return finite_difference_jacobian(
            workspace, psi, fadeev, method, delta,
            settings.finite_difference_scale,
            settings.repair_stationary_points, settings.recovery_sweeps);
    };
    return gss::Newton::solve(
        run_initial_state(workspace, base, settings),
        residual_evaluator,
        jacobian_evaluator,
        direct_linear_config(false),
        newton_config_from(base, damping),
        std::forward<decltype(observer)>(observer));
}

struct ExplicitGResult {
    DenseMatrix G;
    Real max_linear_true_relative_residual{0};
    Index max_linear_iterations{0};
};

/// Form G = A^{-1} B explicitly, column by column.
///
/// This costs n linear solves and n^2 storage and is a verification tool, not
/// a campaign path.  The manuscript already states that explicit formation is
/// unnecessary, since the action of G needs only an application of B followed
/// by a solve with the factorised A; the matrix-free driver is what a campaign
/// should use.
// CHANGED: the same factorisation reuse as run_picard.  n column solves against
// a matrix that was just factorised elsewhere do not need a second factorisation.
[[nodiscard]] inline ExplicitGResult build_explicit_G(
    const CaseWorkspace& workspace,
    gss::PetscLinearSystemController& linear,
    const DenseMatrix& B) {
    const std::size_t n = static_cast<std::size_t>(workspace.poisson.matrix.rows());
    if (B.size() != n * n) throw std::invalid_argument("build_explicit_G B size mismatch");
    if (n > dense_cell_limit) {
        throw std::runtime_error(
            "build_explicit_G refused: forming G explicitly is a small-mesh "
            "verification tool; use the matrix-free operator for a campaign");
    }

    ExplicitGResult result;
    result.G.assign(n * n, Real{0});
    DenseVector rhs(n, Real{0});
    for (std::size_t col = 0; col < n; ++col) {
        for (std::size_t row = 0; row < n; ++row) {
            rhs[row] = B[dense_index(n, row, col)];
        }
        const auto solved = linear.solve(rhs);
        if (!solved.converged()) {
            throw std::runtime_error("A^{-1}B column solve failed: " + solved.reason_string);
        }
        result.max_linear_true_relative_residual = std::max(
            result.max_linear_true_relative_residual,
            static_cast<Real>(solved.true_relative_residual));
        result.max_linear_iterations = std::max(
            result.max_linear_iterations, solved.iterations);
        for (std::size_t row = 0; row < n; ++row) {
            result.G[dense_index(n, row, col)] = solved.solution[row];
        }
    }
    return result;
}

/// Legacy entry point: builds its own controller, as every published campaign did.
[[nodiscard]] inline ExplicitGResult build_explicit_G(
    const CaseWorkspace& workspace,
    const DenseMatrix& B) {
    auto linear = make_picard_controller(workspace);
    return build_explicit_G(workspace, linear, B);
}

[[nodiscard]] inline DenseMatrix relaxed_operator(
    const DenseMatrix& G,
    std::size_t n,
    Real omega) {
    if (G.size() != n * n) throw std::invalid_argument("relaxed_operator size mismatch");
    DenseMatrix result(G.size(), Real{0});
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < n; ++col) {
            const Real identity = row == col ? Real{1} : Real{0};
            result[dense_index(n, row, col)] =
                (Real{1} - omega) * identity + omega * G[dense_index(n, row, col)];
        }
    }
    return result;
}

#ifdef GSMAKER_USE_SLEPC
struct SpectrumSummary {
    Real rho{0};
    std::complex<Real> mu1{};
    std::complex<Real> mu2{};
    PetscInt converged_eigenpairs{0};
    PetscInt eps_iterations{0};
    EPSConvergedReason reason{EPS_CONVERGED_ITERATING};
    Real max_relative_error{0};
    double seconds{0};
};

[[nodiscard]] inline SpectrumSummary spectral_radius(
    const DenseMatrix& matrix,
    std::size_t n,
    const StudySettings& settings) {
    if (matrix.size() != n * n) throw std::invalid_argument("spectral_radius matrix size mismatch");
    auto csr = dense_to_full_csr(matrix, n);
    gss::SlepcEigenConfig config;
    config.eps_type = "krylovschur";
    config.problem_type = EPS_NHEP;
    config.which = EPS_LARGEST_MAGNITUDE;
    config.requested_eigenpairs = std::min<PetscInt>(
        settings.requested_eigenpairs,
        static_cast<PetscInt>(n));
    config.tolerance = settings.eps_tolerance;
    config.max_iterations = settings.eps_max_iterations;
    config.set_from_options = false;
    config.compute_eigenvectors = false;

    gss::SlepcEigenSolverController solver{csr, config};
    const auto result = solver.solve();
    if (result.eigenpairs.empty()) throw std::runtime_error("SLEPc returned no eigenpairs");

    std::vector<gss::SlepcEigenpair> ranked = result.eigenpairs;
    std::sort(ranked.begin(), ranked.end(), [](const auto& lhs, const auto& rhs) {
        return std::abs(lhs.eigenvalue) > std::abs(rhs.eigenvalue);
    });

    SpectrumSummary summary;
    summary.mu1 = {static_cast<Real>(ranked[0].eigenvalue.real()), static_cast<Real>(ranked[0].eigenvalue.imag())};
    summary.rho = static_cast<Real>(std::abs(ranked[0].eigenvalue));
    if (ranked.size() > 1U) {
        summary.mu2 = {static_cast<Real>(ranked[1].eigenvalue.real()), static_cast<Real>(ranked[1].eigenvalue.imag())};
    }
    for (const auto& pair : ranked) {
        summary.max_relative_error = std::max(
            summary.max_relative_error,
            static_cast<Real>(pair.relative_error));
    }
    summary.converged_eigenpairs = result.converged_eigenpairs;
    summary.eps_iterations = result.iterations;
    summary.reason = result.reason;
    summary.seconds = result.solve_seconds;
    return summary;
}

// CHANGED: an overload that reuses a factorisation, for a study that walks a
// trajectory and takes a spectrum at every state.  The default form still
// builds its own, so no existing caller changes.
template<class GradientMethod>
[[nodiscard]] SpectrumSummary spectrum_at_state(
    const CaseWorkspace& workspace,
    gss::PetscLinearSystemController& linear,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real omega = Real{1}) {
    const auto B = dense_B_reference(
        workspace, state, fadeev, method,
        base.solver.source_regularisation_delta,
        settings.finite_difference_scale);
    auto G = build_explicit_G(workspace, linear, B).G;
    if (omega != Real{1}) G = relaxed_operator(G, state.size(), omega);
    return spectral_radius(G, state.size(), settings);
}

template<class GradientMethod>
[[nodiscard]] SpectrumSummary spectrum_at_state(
    const CaseWorkspace& workspace,
    std::span<const Real> state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real omega = Real{1}) {
    auto linear = make_picard_controller(workspace);
    return spectrum_at_state(
        workspace, linear, state, fadeev, method, base, settings, omega);
}
#endif

// The window over which the asymptotic factor is measured.
//
// WHY THIS EXISTS.  A fixed window is wrong in both directions.  Measured over
// 198 relaxed runs with ten ratios, the disagreement between predicted and
// observed factor has median 0.07 per cent on sequences longer than 120
// iterations and median 13.1 per cent on sequences of fifteen or fewer, where
// ten ratios reach back into the transient.  With four ratios the picture
// inverts: the short runs agree to 0.15 per cent while a sequence that ends at
// the rounding floor is measured over ratios that are arithmetic noise.  The
// window therefore scales with the length of the history, and the ratios whose
// residuals lie below the noise floor are discarded rather than averaged.
//
// The spread of the ratios inside the window is reported alongside the factor.
// A sequence that has not settled onto its slowest mode still shows a trend
// across the window, and no choice of window recovers an asymptotic rate from
// it; the caller is told so instead of being given a number that looks like
// one.
struct AsymptoticFactorOptions {
    /// Fraction of the usable history averaged, taken from the end.
    Real window_fraction{Real{1} / Real{3}};
    /// Never average fewer than this many ratios, however short the sequence.
    std::size_t minimum_ratios{3U};
    /// A residual below this multiple of eps times the largest residual is
    /// treated as arithmetic noise and the ratios touching it are discarded.
    Real noise_floor_factor{Real{100}};
    /// Relative spread of the ratios below which the sequence counts as settled.
    Real settled_tolerance{Real{0.05}};

    void validate() const {
        if (!(window_fraction > Real{0}) || window_fraction > Real{1}) {
            throw std::invalid_argument("window_fraction must lie in (0, 1]");
        }
        if (minimum_ratios < 2U) {
            throw std::invalid_argument("minimum_ratios must be at least 2");
        }
        if (!(noise_floor_factor >= Real{0})) {
            throw std::invalid_argument("noise_floor_factor must be non-negative");
        }
        if (!(settled_tolerance > Real{0})) {
            throw std::invalid_argument("settled_tolerance must be positive");
        }
    }
};

struct AsymptoticFactor {
    Real factor{std::numeric_limits<Real>::quiet_NaN()};
    /// Number of ratios that entered the geometric mean.
    std::size_t ratios_used{0U};
    /// Index k of the first ratio used, that is residual[k] / residual[k-1].
    std::size_t first_ratio{0U};
    /// Ratios excluded because a residual lay below the noise floor.
    std::size_t ratios_discarded{0U};
    /// (max - min) / mean over the window.
    Real relative_spread{std::numeric_limits<Real>::quiet_NaN()};
    /// False when the sequence is still approaching its slowest mode.
    bool settled{false};
};

[[nodiscard]] inline AsymptoticFactor asymptotic_factor(
    const std::vector<Real>& residuals,
    const AsymptoticFactorOptions& options = AsymptoticFactorOptions{}) {

    options.validate();

    AsymptoticFactor result;
    if (residuals.size() < 3U) return result;

    const Real largest = *std::max_element(residuals.begin(), residuals.end());
    const Real floor = options.noise_floor_factor *
        std::numeric_limits<Real>::epsilon() * largest;

    std::vector<std::size_t> usable;
    usable.reserve(residuals.size());
    for (std::size_t i = 1U; i < residuals.size(); ++i) {
        if (residuals[i - 1U] > floor && residuals[i] > floor) usable.push_back(i);
    }

    const std::size_t available = residuals.size() - 1U;
    result.ratios_discarded = available - usable.size();
    if (usable.empty()) return result;

    const auto requested = static_cast<std::size_t>(
        std::ceil(options.window_fraction * static_cast<Real>(usable.size())));
    const std::size_t want = std::max(options.minimum_ratios, requested);
    const std::size_t count = std::min(want, usable.size());

    long double log_sum = 0.0L;
    Real smallest_ratio = std::numeric_limits<Real>::max();
    Real largest_ratio = Real{0};
    long double linear_sum = 0.0L;

    for (std::size_t k = usable.size() - count; k < usable.size(); ++k) {
        const std::size_t i = usable[k];
        const Real ratio = residuals[i] / residuals[i - 1U];
        if (!(ratio > Real{0}) || !std::isfinite(static_cast<double>(ratio))) {
            return AsymptoticFactor{};
        }
        log_sum += std::log(static_cast<long double>(ratio));
        linear_sum += static_cast<long double>(ratio);
        smallest_ratio = std::min(smallest_ratio, ratio);
        largest_ratio = std::max(largest_ratio, ratio);
    }

    result.ratios_used = count;
    result.first_ratio = usable[usable.size() - count];
    result.factor = static_cast<Real>(
        std::exp(log_sum / static_cast<long double>(count)));

    const auto mean = static_cast<Real>(
        linear_sum / static_cast<long double>(count));
    result.relative_spread = mean > Real{0}
        ? (largest_ratio - smallest_ratio) / mean
        : std::numeric_limits<Real>::quiet_NaN();
    result.settled = result.relative_spread < options.settled_tolerance;

    return result;
}

/// Fixed-window measurement, retained so that a campaign can reproduce a
/// number published before the proportional window existed.  New work should
/// use asymptotic_factor().
[[nodiscard]] inline Real tail_geometric_factor(
    const std::vector<Real>& residuals,
    std::size_t tail_ratios = 4U,
    std::size_t* ratios_used = nullptr) {
    if (ratios_used != nullptr) *ratios_used = 0U;
    if (residuals.size() < 3U) return std::numeric_limits<Real>::quiet_NaN();
    if (tail_ratios < 1U) tail_ratios = 1U;
    const std::size_t start = residuals.size() > tail_ratios
        ? residuals.size() - tail_ratios
        : 1U;
    long double log_sum = 0.0L;
    std::size_t count = 0U;
    for (std::size_t i = start; i < residuals.size(); ++i) {
        if (residuals[i - 1U] > Real{0} && residuals[i] > Real{0}) {
            const long double ratio = static_cast<long double>(residuals[i] / residuals[i - 1U]);
            log_sum += std::log(ratio);
            ++count;
        }
    }
    if (ratios_used != nullptr) *ratios_used = count;
    return count > 0U
        ? static_cast<Real>(std::exp(log_sum / static_cast<long double>(count)))
        : std::numeric_limits<Real>::quiet_NaN();
}

} // namespace study

// -----------------------------------------------------------------------------
// Additional state-start helpers for the solver-basin and switching studies.
// These remain study-only utilities; no GSMakerLib API is changed.
// -----------------------------------------------------------------------------
namespace study {

// CHANGED: the factorisation is reusable here too, and it matters more in these
// studies than anywhere else: a restart sweep calls this once per sampled state
// and refactorised A every time, so the reported cost of restarting grew with
// the number of restarts for a reason that has nothing to do with the method.
template<class GradientMethod>
[[nodiscard]] inline gss::PicardResult run_picard_from_state(
    const CaseWorkspace& workspace,
    gss::PetscLinearSystemController& linear,
    std::span<const Real> initial_state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real omega,
    auto&& observer) {
    const auto source = make_source_evaluator(workspace, fadeev, method, base, settings);
    return gss::Picard::solve(
        workspace.poisson.matrix,
        workspace.poisson.boundary_rhs,
        initial_state,
        source,
        linear,
        picard_config_from(base, omega),
        std::forward<decltype(observer)>(observer));
}

template<class GradientMethod>
[[nodiscard]] inline gss::PicardResult run_picard_from_state(
    const CaseWorkspace& workspace,
    std::span<const Real> initial_state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    Real omega,
    auto&& observer) {
    auto linear = make_picard_controller(workspace);
    return run_picard_from_state(
        workspace, linear, initial_state, fadeev, method, base,
        StudySettings{}, omega,
        std::forward<decltype(observer)>(observer));
}

// CHANGED: settings-aware overload, so a basin or switching study that enables
// the repair uses the same source map for both linearisations.
template<class GradientMethod>
[[nodiscard]] inline gss::PicardResult run_picard_from_state(
    const CaseWorkspace& workspace,
    std::span<const Real> initial_state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real omega,
    auto&& observer) {
    auto linear = make_picard_controller(workspace);
    return run_picard_from_state(
        workspace, linear, initial_state, fadeev, method, base, settings, omega,
        std::forward<decltype(observer)>(observer));
}

template<class GradientMethod>
[[nodiscard]] inline gss::NewtonResult run_newton_from_state(
    const CaseWorkspace& workspace,
    std::span<const Real> initial_state,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    const gss::Config& base,
    const StudySettings& settings,
    Real damping,
    auto&& observer) {
    const Real delta = base.solver.source_regularisation_delta;
    const auto source = make_source_evaluator(workspace, fadeev, method, base, settings);
    auto residual_evaluator = [&](std::span<const Real> psi) {
        auto value = workspace.poisson.matrix.multiply(psi);
        const auto integrated = source(psi);
        for (std::size_t i = 0; i < value.size(); ++i) {
            value[i] -= integrated[i] + workspace.poisson.boundary_rhs[i];
        }
        return value;
    };
    auto jacobian_evaluator = [&](std::span<const Real> psi) {
        return finite_difference_jacobian(
            workspace, psi, fadeev, method, delta,
            settings.finite_difference_scale,
            settings.repair_stationary_points, settings.recovery_sweeps);
    };
    return gss::Newton::solve(
        initial_state,
        residual_evaluator,
        jacobian_evaluator,
        direct_linear_config(false),
        newton_config_from(base, damping),
        std::forward<decltype(observer)>(observer));
}

[[nodiscard]] inline Real equivalent_picard_source_evaluations(
    const gss::PicardResult& result) noexcept {
    return static_cast<Real>(result.iterations + Index{1});
}

// CHANGED: the cost model follows the assembly.  A coloured Jacobian costs
// 2 * colours source evaluations, not 2n, and the colour count is bounded by
// the maximum distance-2 degree rather than growing with the mesh.  Every
// many-query cost comparison in the campaigns depends on this figure.
[[nodiscard]] inline Real equivalent_newton_source_evaluations(
    const gss::NewtonResult& result,
    std::size_t colours) noexcept {
    return static_cast<Real>(result.iterations + Index{1}) +
           static_cast<Real>(2U * colours) *
               static_cast<Real>(result.jacobian_evaluations);
}

} // namespace study