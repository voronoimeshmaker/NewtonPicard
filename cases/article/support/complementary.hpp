#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <array>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>



#include <petscksp.h>
#include <yaml-cpp/yaml.h>

#include <VoroGradToolkit/Gradients/Common/WLSCore2D.hpp>
#include <VoroGradToolkit/Gradients/GreenGauss/ModifiedGreenGauss2D.hpp>
#include <VoroGradToolkit/Gradients/Polynomial/WLSQR2D.hpp>

#include <adolc/adolc.h>

#include "gsmaker/analysis/spectral/SlepcContext.hpp"
#include "gsmaker/gradient/Gradient.hpp"

#include "common/ConsoleTee.hpp"
#include "common/FadeevNonlinearStudyCommon.hpp"
#include "common/MatrixFreePicardOperator.hpp"

namespace p3official {

namespace ps = study;
namespace sparse = study::sparse;
namespace spectral = study::spectral;
namespace vgt = ::vgt;

using Real = study::Real;
using Index = study::Index;
using DenseVector = std::vector<Real>;

enum class Strategy {
    PicardOnly,
    NewtonOnly,
    MixedRhoInitial,
    MixedRhoEveryIteration
};


struct Settings {
    std::vector<std::string> methods{"WLSQR", "MGG"};
    bool reference_rho{true}, initial_rho{true};
    bool jacobian_timing{false}, floor_activity{false};
    bool bins_by_target{false};
    std::vector<std::string> mesh_families{"cartesian", "hexagonal", "random"};
    Index target_cells{21760};
    Index nominal_nx{256};
    Index nominal_ny{85};
    std::uint32_t mesh_seed{12345U};

    std::size_t sample_count{2000U};
    std::uint64_t sampling_seed{20260914ULL};
    Real alpha_min{Real{1.0e-4}};
    Real alpha_max{Real{9.999e-1}};
    Real theta_min_degrees{Real{-10}};
    Real theta_max_degrees{Real{10}};
    std::array<Real, 6> rho_strata_edges{
        Real{0.70}, Real{0.90}, Real{0.98},
        Real{1.02}, Real{1.10}, Real{1.30}};

    Real calibration_alpha{Real{0.225}};
    Real calibration_theta_degrees{Real{5.7}};

    Real perturbation_amplitude{Real{1.0e-4}};
    Real picard_if_rho_below{Real{1}};
    Real divergence_factor{Real{1.0e8}};
    Real picard_omega{Real{1}};
    Real newton_damping{Real{1}};
    Real reference_relative_tolerance{Real{1.0e-6}};

    // ADOL-C is the project-wide default. Finite differences are available
    // only when source.derivative explicitly requests them. There is no
    // automatic fallback from ADOL-C.
    std::string source_derivative{"adolc"};

    std::size_t max_iterations{5000U};
    Real tolerance_relative{Real{1.0e-10}};
    Real tolerance_absolute{Real{1.0e-12}};

    PetscInt requested_eigenpairs{8};
    PetscReal eps_tolerance{1.0e-11};
    PetscInt eps_max_iterations{3000};

    PetscReal krylov_relative_tolerance{1.0e-2};
    PetscInt krylov_max_iterations{50};

    bool compute_unexpected_picard_final_rho{true};
    std::size_t progress_every{100U};
};

[[nodiscard]] inline Settings load_settings(const std::string& path) {
    Settings s;
    const YAML::Node root = YAML::LoadFile(path);
    if(root["methods"]) s.methods=root["methods"].as<std::vector<std::string>>();
    if(const auto d=root["diagnostics"]) {
        s.reference_rho=d["reference_rho"].as<bool>(true);
        s.initial_rho=d["initial_rho"].as<bool>(true);
        s.jacobian_timing=d["jacobian_assembly"].as<bool>(false);
        s.floor_activity=d["regularisation_floor"].as<bool>(false);
    }
    if(root["analysis"]) s.bins_by_target=root["analysis"]["rho_bins_by_target"].as<bool>(false);
    if(s.methods.empty()) throw std::invalid_argument("methods must not be empty");
    for(const auto& m:s.methods) if(m!="WLSQR" && m!="MGG") throw std::invalid_argument("unknown gradient method");

    if (const auto mesh = root["mesh"]) {
        if (mesh["families"])
            s.mesh_families = mesh["families"].as<std::vector<std::string>>();
        if (mesh["target_cells"])
            s.target_cells = static_cast<Index>(mesh["target_cells"].as<long long>());
        if (mesh["nominal_nx"])
            s.nominal_nx = static_cast<Index>(mesh["nominal_nx"].as<long long>());
        if (mesh["nominal_ny"])
            s.nominal_ny = static_cast<Index>(mesh["nominal_ny"].as<long long>());
        if (mesh["random_seed"])
            s.mesh_seed = mesh["random_seed"].as<std::uint32_t>();
    }

    if(const auto mesh=root["mesh"]) {
        if(mesh["n_per_direction"]) s.nominal_nx=mesh["n_per_direction"].as<Index>();
        if(mesh["n_per_direction_y"]) s.nominal_ny=mesh["n_per_direction_y"].as<Index>();
        if(s.nominal_nx<=0 || s.nominal_ny<=0 || s.target_cells!=s.nominal_nx*s.nominal_ny ||
           (mesh["nominal_nx"] && mesh["nominal_nx"].as<Index>()!=s.nominal_nx) ||
           (mesh["nominal_ny"] && mesh["nominal_ny"].as<Index>()!=s.nominal_ny) ||
           (mesh["n_generators"] && mesh["n_generators"].as<Index>()!=s.target_cells))
            throw std::invalid_argument("inconsistent mesh: target_cells must equal nx*ny; nominal sizes and n_generators must agree");
    }
    if (const auto sampling = root["sampling"]) {
        if (sampling["sample_count"])
            s.sample_count = sampling["sample_count"].as<std::size_t>();
        if (sampling["seed"])
            s.sampling_seed = sampling["seed"].as<std::uint64_t>();
        if (sampling["alpha_min"])
            s.alpha_min = static_cast<Real>(sampling["alpha_min"].as<double>());
        if (sampling["alpha_max"])
            s.alpha_max = static_cast<Real>(sampling["alpha_max"].as<double>());
        if (sampling["theta_min_degrees"])
            s.theta_min_degrees =
                static_cast<Real>(sampling["theta_min_degrees"].as<double>());
        if (sampling["theta_max_degrees"])
            s.theta_max_degrees =
                static_cast<Real>(sampling["theta_max_degrees"].as<double>());
        if (sampling["rho_strata_edges"]) {
            const auto edges = sampling["rho_strata_edges"].as<std::vector<double>>();
            if (edges.size() != s.rho_strata_edges.size())
                throw std::invalid_argument(
                    "sampling.rho_strata_edges must contain exactly 6 values");
            for (std::size_t i = 0; i < edges.size(); ++i)
                s.rho_strata_edges[i] = static_cast<Real>(edges[i]);
        }
    }

    if (const auto calibration = root["calibration"]) {
        if (calibration["alpha"])
            s.calibration_alpha =
                static_cast<Real>(calibration["alpha"].as<double>());
        if (calibration["theta_degrees"])
            s.calibration_theta_degrees =
                static_cast<Real>(calibration["theta_degrees"].as<double>());
    }

    if (const auto source = root["source"]) {
        if (source["derivative"])
            s.source_derivative = source["derivative"].as<std::string>();
    }

    if (const auto solver = root["solver"]) {
        if (solver["max_iter"])
            s.max_iterations = solver["max_iter"].as<std::size_t>();
        if (solver["tolerance_relative"])
            s.tolerance_relative =
                static_cast<Real>(solver["tolerance_relative"].as<double>());
        if (solver["tolerance_absolute"])
            s.tolerance_absolute =
                static_cast<Real>(solver["tolerance_absolute"].as<double>());
        if (solver["relaxation_omega"])
            s.picard_omega =
                static_cast<Real>(solver["relaxation_omega"].as<double>());
        if (solver["perturbation_amplitude"])
            s.perturbation_amplitude =
                static_cast<Real>(solver["perturbation_amplitude"].as<double>());
    }

    if (const auto spectral_node = root["spectral"]) {
        if (spectral_node["requested_eigenpairs"])
            s.requested_eigenpairs =
                static_cast<PetscInt>(spectral_node["requested_eigenpairs"].as<int>());
        if (spectral_node["tolerance"])
            s.eps_tolerance =
                static_cast<PetscReal>(spectral_node["tolerance"].as<double>());
        if (spectral_node["max_iterations"])
            s.eps_max_iterations =
                static_cast<PetscInt>(spectral_node["max_iterations"].as<int>());
    }

    if (const auto study_node = root["study"]) {
        if (study_node["newton_damping"])
            s.newton_damping =
                static_cast<Real>(study_node["newton_damping"].as<double>());
        if (study_node["reference_relative_tolerance"])
            s.reference_relative_tolerance =
                static_cast<Real>(study_node["reference_relative_tolerance"].as<double>());
        if (study_node["progress_every"])
            s.progress_every = study_node["progress_every"].as<std::size_t>();
    }

    if (const auto comparison = root["comparison"]) {
        if (comparison["picard_if_rho_below"])
            s.picard_if_rho_below =
                static_cast<Real>(comparison["picard_if_rho_below"].as<double>());
        if (comparison["divergence_factor"])
            s.divergence_factor =
                static_cast<Real>(comparison["divergence_factor"].as<double>());
        if (comparison["compute_unexpected_picard_final_rho"])
            s.compute_unexpected_picard_final_rho =
                comparison["compute_unexpected_picard_final_rho"].as<bool>();
    }

    if (const auto krylov = root["krylov"]) {
        if (krylov["relative_tolerance"])
            s.krylov_relative_tolerance =
                static_cast<PetscReal>(krylov["relative_tolerance"].as<double>());
        if (krylov["max_iterations"])
            s.krylov_max_iterations =
                static_cast<PetscInt>(krylov["max_iterations"].as<int>());
    }

    if (s.mesh_families.empty())
        throw std::invalid_argument("mesh.families must not be empty");
    if (s.sample_count == 0U)
        throw std::invalid_argument("sampling.sample_count must be positive");
    if (!(s.alpha_min > Real{0} && s.alpha_min < s.alpha_max && s.alpha_max < Real{1}))
        throw std::invalid_argument("sampling alpha range must satisfy 0 < alpha_min < alpha_max < 1");
    if (!(s.theta_min_degrees < s.theta_max_degrees))
        throw std::invalid_argument("invalid theta range");
    for (std::size_t i = 0; i + 1U < s.rho_strata_edges.size(); ++i) {
        if (!(s.rho_strata_edges[i] > Real{0} &&
              s.rho_strata_edges[i] < s.rho_strata_edges[i + 1U]))
            throw std::invalid_argument(
                "sampling.rho_strata_edges must be positive and strictly increasing");
    }
    if (!(s.picard_if_rho_below > Real{0}))
        throw std::invalid_argument("comparison.picard_if_rho_below must be positive");
    if (!(s.reference_relative_tolerance > Real{0}))
        throw std::invalid_argument("study.reference_relative_tolerance must be positive");
    if (!(s.perturbation_amplitude > Real{0}))
        throw std::invalid_argument("solver.perturbation_amplitude must be positive");
    if (s.max_iterations == 0U)
        throw std::invalid_argument("solver.max_iter must be positive");
    if (s.source_derivative != "adolc" &&
        s.source_derivative != "finite_difference")
        throw std::invalid_argument(
            "source.derivative must be 'adolc' or the explicitly requested 'finite_difference'");

    return s;
}

[[nodiscard]] inline const char* strategy_name(Strategy strategy) noexcept {
    switch (strategy) {
        case Strategy::PicardOnly: return "picard_only";
        case Strategy::NewtonOnly: return "newton_only";
        case Strategy::MixedRhoInitial: return "mixed_rho_initial";
        case Strategy::MixedRhoEveryIteration: return "mixed_rho_every_iteration";
    }
    return "unknown";
}

// Compatibility wrapper for the current and immediately preceding Common APIs.
template<class GradientMethod>
[[nodiscard]] DenseVector evaluate_base_source(
    const ps::CaseWorkspace& workspace,
    std::span<const Real> psi,
    const gss::Fadeev& fadeev,
    const GradientMethod& method,
    Real delta,
    const ps::StudySettings& study_settings) {

    if constexpr (requires {
        ps::evaluate_integrated_source(
            workspace,
            psi,
            fadeev,
            method,
            delta,
            study_settings.repair_stationary_points,
            study_settings.recovery_sweeps);
    }) {
        return ps::evaluate_integrated_source(
            workspace,
            psi,
            fadeev,
            method,
            delta,
            study_settings.repair_stationary_points,
            study_settings.recovery_sweeps);
    } else {
        return ps::evaluate_integrated_source(
            workspace, psi, fadeev, method, delta);
    }
}

[[nodiscard]] inline DenseVector make_compensator(
    const ps::CaseWorkspace& workspace,
    std::span<const Real> reference_state,
    std::span<const Real> base_source_at_reference,
    Real scale) {

    if (reference_state.size() != base_source_at_reference.size())
        throw std::invalid_argument("make_compensator: size mismatch");

    const auto Apsi = workspace.poisson.matrix.multiply(reference_state);
    const auto& boundary = workspace.poisson.boundary_rhs;

    if (Apsi.size() != boundary.size())
        throw std::runtime_error(
            "make_compensator: Poisson matrix/boundary size mismatch");

    DenseVector q(reference_state.size(), Real{0});
    for (std::size_t i = 0; i < q.size(); ++i)
        q[i] = Apsi[i] - boundary[i] - scale * base_source_at_reference[i];

    return q;
}

template<class GradientMethod>
class ModifiedSource final {
public:
    ModifiedSource(
        const ps::CaseWorkspace& workspace,
        const gss::Fadeev& fadeev,
        const GradientMethod& method,
        Real delta,
        const ps::StudySettings& study_settings,
        Real scale,
        const DenseVector& compensator, bool observe_floor=false)
        : workspace_(workspace),
          fadeev_(fadeev),
          method_(method),
          delta_(delta),
          study_settings_(study_settings),
          scale_(scale),
          compensator_(compensator), observe_floor_(observe_floor) {}

    [[nodiscard]] DenseVector operator()(std::span<const Real> psi) const {
        auto source = observe_floor_
            ? ps::evaluate_integrated_source(workspace_, psi, fadeev_, method_, delta_,
                study_settings_.repair_stationary_points, study_settings_.recovery_sweeps, &last_diagnostics_)
            : evaluate_base_source(workspace_, psi, fadeev_, method_, delta_, study_settings_);

        if (source.size() != compensator_.size())
            throw std::runtime_error("ModifiedSource: source/compensator size mismatch");

        for (std::size_t i = 0; i < source.size(); ++i)
            source[i] = scale_ * source[i] + compensator_[i];

        return source;
    }

    bool floor_active() const noexcept { return observe_floor_ && last_diagnostics_.regularised_cells>0; }
    [[nodiscard]] const gss::Fadeev& fadeev() const noexcept {
        return fadeev_;
    }

    [[nodiscard]] Real scale() const noexcept {
        return scale_;
    }

private:
    const ps::CaseWorkspace& workspace_;
    const gss::Fadeev& fadeev_;
    const GradientMethod& method_;
    Real delta_;
    const ps::StudySettings& study_settings_;
    Real scale_;
    const DenseVector& compensator_;
    bool observe_floor_{false};
    mutable ps::SourceDiagnostics last_diagnostics_;
};

// -----------------------------------------------------------------------------
// Sparse storage and ADOL-C source Jacobian
// -----------------------------------------------------------------------------

[[nodiscard]] inline std::vector<Index> csr_row_columns(
    const gss::CsrMatrix& matrix,
    std::size_t row) {

    const auto starts = matrix.row_start();
    const auto columns = matrix.column();
    const auto begin = static_cast<std::size_t>(starts[row]);
    const auto end = static_cast<std::size_t>(starts[row + 1U]);
    return std::vector<Index>(columns.begin() + begin, columns.begin() + end);
}

inline void require_matrix_ring_pattern(
    const gss::CsrMatrix& matrix,
    const gss::MeshTopology& topology,
    std::size_t radius,
    std::string_view name) {

    const auto expected = sparse::closed_rings(topology, radius);
    if (static_cast<std::size_t>(matrix.rows()) != expected.size())
        throw std::runtime_error(
            std::string{name} + ": row count does not match mesh adjacency");

    for (std::size_t row = 0; row < expected.size(); ++row) {
        const auto actual = csr_row_columns(matrix, row);
        if (actual != expected[row]) {
            throw std::runtime_error(
                std::string{name} +
                ": CSR allocation does not match the requested mesh adjacency ring at row " +
                std::to_string(row));
        }
    }
}

inline void require_workspace_sparse_storage(
    const ps::CaseWorkspace& workspace) {

    const auto& topology = workspace.mesh.topology();

    // The Poisson operator is a two-point finite-volume operator, hence its
    // matrix must be exactly diagonal plus direct face neighbours.
    require_matrix_ring_pattern(
        workspace.poisson.matrix, topology, 1U, "Poisson A");

    // B is allocated from the configured graph ring. With the production
    // setting radius=1 this is the same face-neighbour adjacency matrix.
    const auto expected_b =
        sparse::closed_rings(topology, workspace.jacobian_pattern.radius);

    if (workspace.jacobian_pattern.row_start.size() != expected_b.size() + 1U)
        throw std::runtime_error(
            "source-Jacobian CSR row allocation does not match the mesh");

    for (std::size_t row = 0; row < expected_b.size(); ++row) {
        const auto begin = static_cast<std::size_t>(
            workspace.jacobian_pattern.row_start[row]);
        const auto end = static_cast<std::size_t>(
            workspace.jacobian_pattern.row_start[row + 1U]);
        const std::vector<Index> actual(
            workspace.jacobian_pattern.column.begin() + begin,
            workspace.jacobian_pattern.column.begin() + end);
        if (actual != expected_b[row]) {
            throw std::runtime_error(
                "source-Jacobian CSR allocation does not match the configured mesh ring at row " +
                std::to_string(row));
        }
    }
}

[[nodiscard]] inline int checked_adolc_dimension(
    std::size_t value,
    std::string_view what) {

    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(
            std::string{"ADOL-C dimension exceeds int range: "} +
            std::string{what});
    return static_cast<int>(value);
}

template<class GradientMethod>
class SourceJacobianAssembler final {
public:
    SourceJacobianAssembler(
        const ps::CaseWorkspace& workspace,
        const GradientMethod& method,
        Real regularisation_delta,
        const ps::StudySettings& study_settings,
        std::string derivative_method)
        : pattern_(workspace.jacobian_pattern),
          volumes_(workspace.volumes),
          regularisation_delta_(regularisation_delta),
          finite_difference_scale_(study_settings.finite_difference_scale),
          derivative_method_(std::move(derivative_method)) {

        require_workspace_sparse_storage(workspace);

        if (derivative_method_ == "finite_difference") {
            // This path exists only for a formally requested comparison.
            return;
        }
        if (derivative_method_ != "adolc")
            throw std::invalid_argument(
                "SourceJacobianAssembler: unsupported derivative method");

        if (study_settings.repair_stationary_points)
            throw std::invalid_argument(
                "ADOL-C source Jacobian does not silently differentiate the "
                "non-smooth stationary-point repair; disable repair or implement "
                "an explicit active-compatible rule");

        if (!(regularisation_delta_ > Real{0}) ||
            !std::isfinite(static_cast<double>(regularisation_delta_)))
            throw std::invalid_argument(
                "ADOL-C source Jacobian requires a finite positive regularisation delta");

        build_gradient_operator(workspace, method);
        verify_gradient_operator(workspace, method);
        build_seed_matrix();
        trace_base_source(workspace.exact);
    }

    SourceJacobianAssembler(const SourceJacobianAssembler&) = delete;
    SourceJacobianAssembler& operator=(const SourceJacobianAssembler&) = delete;

    template<class Source>
    [[nodiscard]] std::vector<Real> values(
        std::span<const Real> state,
        const Source& source) {

        if (state.size() != pattern_.cells)
            throw std::invalid_argument(
                "SourceJacobianAssembler: state size does not match mesh");

        if (derivative_method_ == "finite_difference") {
            return sparse::assemble_values(
                pattern_, state, source, finite_difference_scale_);
        }

        return adolc_values(state, source.fadeev(), source.scale());
    }

    [[nodiscard]] const std::string& derivative_method() const noexcept {
        return derivative_method_;
    }

    [[nodiscard]] Real gradient_operator_relative_defect() const noexcept {
        return gradient_operator_relative_defect_;
    }

private:
    void apply_gradient_operator(
        std::span<const Real> field,
        DenseVector& gx,
        DenseVector& gy) const {

        gx.assign(pattern_.cells, Real{0});
        gy.assign(pattern_.cells, Real{0});

        for (std::size_t row = 0; row < pattern_.cells; ++row) {
            const auto begin = static_cast<std::size_t>(pattern_.row_start[row]);
            const auto end = static_cast<std::size_t>(pattern_.row_start[row + 1U]);
            for (std::size_t k = begin; k < end; ++k) {
                const auto col = static_cast<std::size_t>(pattern_.column[k]);
                gx[row] += gradient_x_[k] * field[col];
                gy[row] += gradient_y_[k] * field[col];
            }
        }
    }

    template<class GradientMethod2>
    void build_gradient_operator(
        const ps::CaseWorkspace& workspace,
        const GradientMethod2& method) {

        gradient_x_.assign(pattern_.nnz(), Real{0});
        gradient_y_.assign(pattern_.nnz(), Real{0});

        DenseVector probe(pattern_.cells, Real{0});

        // The reconstruction is linear for fixed geometry. A graph-coloured
        // exact superposition recovers its local coefficients without any
        // finite-difference quotient. The conflict colouring already owned by
        // the source-Jacobian pattern guarantees at most one active column of
        // a colour in each row.
        for (std::size_t colour = 0; colour < pattern_.colours; ++colour) {
            for (const Index member : pattern_.colour_members[colour])
                probe[static_cast<std::size_t>(member)] = Real{1};

            const auto reconstructed =
                gss::Gradient::compute(workspace.gradients, probe, method);

            if (reconstructed.diagnostics.failed_volumes != 0U)
                throw std::runtime_error(
                    "gradient-operator extraction failed in VoroGradToolkit");

            for (std::size_t row = 0; row < pattern_.cells; ++row) {
                const auto begin = static_cast<std::size_t>(pattern_.row_start[row]);
                const auto end = static_cast<std::size_t>(pattern_.row_start[row + 1U]);
                for (std::size_t k = begin; k < end; ++k) {
                    const auto col = static_cast<std::size_t>(pattern_.column[k]);
                    if (static_cast<std::size_t>(pattern_.colour[col]) != colour)
                        continue;
                    gradient_x_[k] = reconstructed.field.x[row];
                    gradient_y_[k] = reconstructed.field.y[row];
                }
            }

            for (const Index member : pattern_.colour_members[colour])
                probe[static_cast<std::size_t>(member)] = Real{0};
        }
    }

    template<class GradientMethod2>
    void verify_gradient_operator(
        const ps::CaseWorkspace& workspace,
        const GradientMethod2& method) {

        DenseVector field(pattern_.cells, Real{0});
        for (std::size_t i = 0; i < field.size(); ++i) {
            const double q = static_cast<double>(i + 1U);
            field[i] = static_cast<Real>(
                std::sin(0.173 * q) + 0.37 * std::cos(0.419 * q));
        }

        const auto native =
            gss::Gradient::compute(workspace.gradients, field, method);

        DenseVector gx, gy;
        apply_gradient_operator(field, gx, gy);

        long double numerator = 0.0L;
        long double denominator = 0.0L;
        for (std::size_t i = 0; i < field.size(); ++i) {
            const long double dx =
                static_cast<long double>(gx[i]) -
                static_cast<long double>(native.field.x[i]);
            const long double dy =
                static_cast<long double>(gy[i]) -
                static_cast<long double>(native.field.y[i]);
            numerator += dx * dx + dy * dy;
            denominator +=
                static_cast<long double>(native.field.x[i]) *
                static_cast<long double>(native.field.x[i]) +
                static_cast<long double>(native.field.y[i]) *
                static_cast<long double>(native.field.y[i]);
        }

        gradient_operator_relative_defect_ = static_cast<Real>(
            std::sqrt(numerator / std::max(denominator, 1.0e-300L)));

        if (!std::isfinite(static_cast<double>(gradient_operator_relative_defect_)) ||
            gradient_operator_relative_defect_ > Real{1.0e-10})
            throw std::runtime_error(
                "ADOL-C gradient operator does not reproduce the native "
                "VoroGradToolkit reconstruction on the configured adjacency stencil");
    }

    void build_seed_matrix() {
        const std::size_t directions = pattern_.colours;
        const std::size_t independent_variables = pattern_.cells + 1U;
        if (directions == 0U)
            throw std::runtime_error(
                "ADOL-C source Jacobian received an empty graph colouring");

        // The extra independent variable is the regularisation threshold.
        // Its seed is identically zero, so AD differentiates the source with
        // the numerical regularisation scale frozen at the current state.
        seed_data_.assign(independent_variables * directions, 0.0);
        seed_rows_.resize(independent_variables);
        derivative_data_.assign(pattern_.cells * directions, 0.0);
        derivative_rows_.resize(pattern_.cells);

        for (std::size_t i = 0; i < independent_variables; ++i)
            seed_rows_[i] = seed_data_.data() + i * directions;

        for (std::size_t i = 0; i < pattern_.cells; ++i) {
            derivative_rows_[i] = derivative_data_.data() + i * directions;
            const auto colour = static_cast<std::size_t>(pattern_.colour[i]);
            seed_rows_[i][colour] = 1.0;
        }
    }

    [[nodiscard]] Real frozen_regularisation_threshold(
        std::span<const Real> state) const {

        DenseVector gx, gy;
        apply_gradient_operator(state, gx, gy);

        Real largest = Real{0};
        for (std::size_t i = 0; i < state.size(); ++i) {
            const Real g2 = gx[i] * gx[i] + gy[i] * gy[i];
            largest = std::max(largest, g2);
        }

        const Real threshold =
            regularisation_delta_ * regularisation_delta_ * largest;

        if (!(threshold > Real{0}) ||
            !std::isfinite(static_cast<double>(threshold)))
            throw std::runtime_error(
                "ADOL-C source Jacobian cannot linearise a zero/non-finite gradient scale");

        return threshold;
    }

    void trace_base_source(std::span<const Real> reference_state) {
        const std::size_t n = pattern_.cells;
        if (reference_state.size() != n)
            throw std::invalid_argument(
                "ADOL-C trace reference state does not match the mesh");

        const Real threshold = frozen_regularisation_threshold(reference_state);

        trace_on(tape_tag_);

        // Psi occupies [0,n); the last independent is the threshold.  The
        // threshold is evaluated at every requested state but receives zero
        // directional seeds.  This keeps the regularisation safeguard local
        // to the linearisation and prevents its mesh-wide max() scale from
        // creating artificial global Jacobian storage.
        std::vector<adouble> x(n + 1U);
        std::vector<adouble> free_function(n);
        std::vector<adouble> integrated_source(n);
        std::vector<double> traced_output(n, 0.0);

        for (std::size_t i = 0; i < n; ++i)
            x[i] <<= static_cast<double>(reference_state[i]);
        x[n] <<= static_cast<double>(threshold);

        for (std::size_t i = 0; i < n; ++i) {
            adouble exponent = -2.0 * x[i];
            adouble cap = 700.0;
            adouble capped;
            condassign(capped, exponent - cap, cap, exponent);
            using std::exp;
            free_function[i] = -0.5 * exp(capped);
        }

        for (std::size_t row = 0; row < n; ++row) {
            adouble gx_psi = 0.0;
            adouble gy_psi = 0.0;
            adouble gx_f = 0.0;
            adouble gy_f = 0.0;

            const auto begin = static_cast<std::size_t>(pattern_.row_start[row]);
            const auto end = static_cast<std::size_t>(pattern_.row_start[row + 1U]);
            for (std::size_t k = begin; k < end; ++k) {
                const auto col = static_cast<std::size_t>(pattern_.column[k]);
                const double cx = static_cast<double>(gradient_x_[k]);
                const double cy = static_cast<double>(gradient_y_[k]);
                gx_psi += cx * x[col];
                gy_psi += cy * x[col];
                gx_f += cx * free_function[col];
                gy_f += cy * free_function[col];
            }

            const adouble numerator = gx_f * gx_psi + gy_f * gy_psi;
            const adouble denominator =
                gx_psi * gx_psi + gy_psi * gy_psi + x[n];

            integrated_source[row] =
                -static_cast<double>(volumes_[row]) * numerator / denominator;
        }

        for (std::size_t i = 0; i < n; ++i)
            integrated_source[i] >>= traced_output[i];

        trace_off();
    }

    [[nodiscard]] std::vector<Real> adolc_values(
        std::span<const Real> state,
        const gss::Fadeev& fadeev,
        Real outer_scale) {

        for (const Real value : state) {
            if (!std::isfinite(static_cast<double>(value)))
                throw std::invalid_argument(
                    "ADOL-C source Jacobian received a non-finite state");
        }

        const Real threshold = frozen_regularisation_threshold(state);
        const std::size_t n = pattern_.cells;
        const std::size_t independent_variables = n + 1U;

        auto& point = point_buffer_;
        auto& output = output_buffer_;
        point.resize(independent_variables);
        output.resize(n);
        for (std::size_t i = 0; i < n; ++i)
            point[i] = static_cast<double>(state[i]);
        point[n] = static_cast<double>(threshold);

        std::fill(derivative_data_.begin(), derivative_data_.end(), 0.0);

        const int rc = fov_forward(
            tape_tag_,
            checked_adolc_dimension(n, "outputs"),
            checked_adolc_dimension(independent_variables, "inputs"),
            checked_adolc_dimension(pattern_.colours, "directions"),
            point.data(),
            seed_rows_.data(),
            output.data(),
            derivative_rows_.data());

        if (rc < 0)
            throw std::runtime_error(
                "ADOL-C fov_forward failed while assembling the source Jacobian");

        // F is linear in Fadeev's nonlinearity scale and ModifiedSource is
        // linear in its outer manufactured-family scale.  The compensator and
        // manufactured q are fixed data, so neither contributes to B.
        const Real derivative_scale =
            outer_scale * fadeev.nonlinearity_scale();

        std::vector<Real> values(pattern_.nnz(), Real{0});
        for (std::size_t row = 0; row < n; ++row) {
            const auto begin = static_cast<std::size_t>(pattern_.row_start[row]);
            const auto end = static_cast<std::size_t>(pattern_.row_start[row + 1U]);
            for (std::size_t k = begin; k < end; ++k) {
                const auto col = static_cast<std::size_t>(pattern_.column[k]);
                const auto colour = static_cast<std::size_t>(pattern_.colour[col]);
                const double value = derivative_rows_[row][colour];
                if (!std::isfinite(value))
                    throw std::runtime_error(
                        "ADOL-C produced a non-finite source-Jacobian entry");
                values[k] = derivative_scale * static_cast<Real>(value);
            }
        }

        return values;
    }

    const sparse::Pattern& pattern_;
    DenseVector volumes_;
    Real regularisation_delta_{Real{0}};
    Real finite_difference_scale_{Real{1}};
    std::string derivative_method_;

    std::vector<Real> gradient_x_;
    std::vector<Real> gradient_y_;
    Real gradient_operator_relative_defect_{std::numeric_limits<Real>::quiet_NaN()};

    std::vector<double> point_buffer_, output_buffer_;
    std::vector<double> seed_data_;
    std::vector<double*> seed_rows_;
    std::vector<double> derivative_data_;
    std::vector<double*> derivative_rows_;

    // The 10,000-case programs run serially. One tape is built per invariant
    // mesh/method Jacobian assembler and is consumed before another assembler
    // can reuse the tag.
    static constexpr short tape_tag_{23117};
};

[[nodiscard]] inline DenseVector residual_from_source(
    const ps::CaseWorkspace& workspace,
    std::span<const Real> nonlinear,
    std::span<const Real> state) {

    auto result = workspace.poisson.matrix.multiply(state);

    if (result.size() != nonlinear.size() ||
        result.size() != workspace.poisson.boundary_rhs.size())
        throw std::runtime_error("residual: size mismatch");

    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] -= nonlinear[i] + workspace.poisson.boundary_rhs[i];

    return result;
}

[[nodiscard]] inline DenseVector perturbed_state(
    const ps::CaseWorkspace& workspace,
    std::span<const Real> reference_state,
    Real amplitude) {

    DenseVector state(reference_state.begin(), reference_state.end());
    const auto& topology = workspace.mesh.topology();

    DenseVector shape(state.size(), Real{0});
    Real maximum = Real{0};

    for (std::size_t i = 0; i < state.size(); ++i) {
        const double x = static_cast<double>(topology.coord_x[i]);
        const double y = static_cast<double>(topology.coord_y[i]);
        const Real value = static_cast<Real>(
            std::sin(0.731 * x + 1.173 * y) +
            0.5 * std::cos(1.337 * x - 0.619 * y));
        shape[i] = value;
        maximum = std::max(
            maximum,
            static_cast<Real>(std::abs(static_cast<double>(value))));
    }

    if (!(maximum > Real{0}))
        throw std::runtime_error("perturbed_state: degenerate perturbation shape");

    for (std::size_t i = 0; i < state.size(); ++i)
        state[i] += amplitude * shape[i] / maximum;

    return state;
}

class PreconditionedNewtonSolver final {
public:
    PreconditionedNewtonSolver(
        const gss::CsrMatrix& a,
        const gss::CsrMatrix& j_initial,
        const sparse::Pattern& b_pattern,
        bool a_is_spd,
        PetscReal relative_tolerance,
        PetscInt max_iterations) {

        const auto rows = static_cast<PetscInt>(a.rows());
        if (rows <= 0)
            throw std::invalid_argument("PreconditionedNewtonSolver: empty matrix");
        if (j_initial.rows() != a.rows())
            throw std::invalid_argument(
                "PreconditionedNewtonSolver: A/J size mismatch");

        cells_ = static_cast<std::size_t>(rows);

        copy_csr(a, a_row_, a_col_, a_val_);
        copy_csr(j_initial, j_row_, j_col_, j_val_);
        j_base_.assign(j_val_.size(), PetscScalar{0});
        b_to_j_.resize(b_pattern.nnz());
        for (std::size_t row = 0; row < cells_; ++row) {
            const auto begin = j_col_.begin() + j_row_[row];
            const auto end = j_col_.begin() + j_row_[row + 1U];
            const auto position = [&](Index col) {
                const auto it = std::lower_bound(begin, end, static_cast<PetscInt>(col));
                if (it == end || *it != static_cast<PetscInt>(col))
                    throw std::runtime_error("J pattern does not contain A/B entry");
                return static_cast<std::size_t>(it - j_col_.begin());
            };
            for (auto k = a_row_[row]; k < a_row_[row + 1U]; ++k)
                j_base_[position(a_col_[k])] = a_val_[k];
            for (auto k = b_pattern.row_start[row]; k < b_pattern.row_start[row + 1U]; ++k)
                b_to_j_[k] = position(b_pattern.column[k]);
        }

        check(MatCreateSeqAIJWithArrays(
                  PETSC_COMM_SELF, rows, rows,
                  a_row_.data(), a_col_.data(), a_val_.data(), &a_matrix_),
              "creating A");

        check(MatCreateSeqAIJWithArrays(
                  PETSC_COMM_SELF, rows, rows,
                  j_row_.data(), j_col_.data(), j_val_.data(), &j_matrix_),
              "creating J");

        const auto start = std::chrono::steady_clock::now();

        check(KSPCreate(PETSC_COMM_SELF, &ksp_), "creating KSP");
        check(KSPSetOperators(ksp_, j_matrix_, a_matrix_), "setting operators");
        check(KSPSetType(ksp_, KSPGMRES), "selecting GMRES");
        check(KSPSetTolerances(
                  ksp_,
                  relative_tolerance,
                  PETSC_DEFAULT,
                  PETSC_DEFAULT,
                  max_iterations),
              "setting KSP tolerances");
        check(KSPSetErrorIfNotConverged(ksp_, PETSC_FALSE),
              "allowing non-fatal KSP non-convergence");

        PC pc = nullptr;
        check(KSPGetPC(ksp_, &pc), "retrieving PC");
        check(PCSetType(pc, a_is_spd ? PCCHOLESKY : PCLU), "selecting PC");
        check(KSPSetReusePreconditioner(ksp_, PETSC_TRUE), "reusing PC");
        check(KSPSetUp(ksp_), "factorising A");

        factorisation_seconds_ =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();

        check(MatCreateVecs(j_matrix_, &solution_, &rhs_), "creating vectors");
    }

    PreconditionedNewtonSolver(const PreconditionedNewtonSolver&) = delete;
    PreconditionedNewtonSolver& operator=(const PreconditionedNewtonSolver&) = delete;
    PreconditionedNewtonSolver(PreconditionedNewtonSolver&&) = delete;
    PreconditionedNewtonSolver& operator=(PreconditionedNewtonSolver&&) = delete;

    ~PreconditionedNewtonSolver() {
        if (rhs_ != nullptr) VecDestroy(&rhs_);
        if (solution_ != nullptr) VecDestroy(&solution_);
        if (ksp_ != nullptr) KSPDestroy(&ksp_);
        if (j_matrix_ != nullptr) MatDestroy(&j_matrix_);
        if (a_matrix_ != nullptr) MatDestroy(&a_matrix_);
    }

    void update_source_jacobian_values(std::span<const Real> b_values) {
        if (b_values.size() != b_to_j_.size())
            throw std::invalid_argument("source Jacobian size changed");
        std::copy(j_base_.begin(), j_base_.end(), j_val_.begin());
        for (std::size_t k = 0; k < b_values.size(); ++k)
            j_val_[b_to_j_[k]] -= static_cast<PetscScalar>(b_values[k]);
        check(MatAssemblyBegin(j_matrix_, MAT_FINAL_ASSEMBLY), "assembling J begin");
        check(MatAssemblyEnd(j_matrix_, MAT_FINAL_ASSEMBLY), "assembling J end");
    }

    [[nodiscard]] Index solve(
        std::span<const Real> right_hand_side,
        DenseVector& out) {

        if (right_hand_side.size() != cells_)
            throw std::invalid_argument(
                "PreconditionedNewtonSolver: RHS size mismatch");

        PetscScalar* rhs_values = nullptr;
        check(VecGetArray(rhs_, &rhs_values), "accessing RHS");
        for (std::size_t i = 0; i < cells_; ++i)
            rhs_values[i] = static_cast<PetscScalar>(right_hand_side[i]);
        check(VecRestoreArray(rhs_, &rhs_values), "releasing RHS");

        check(KSPSolve(ksp_, rhs_, solution_), "solving Newton correction");
        check(KSPGetConvergedReason(ksp_, &last_reason_), "reading KSP reason");

        PetscInt iterations = 0;
        check(KSPGetIterationNumber(ksp_, &iterations), "reading KSP iterations");

        const PetscScalar* solution_values = nullptr;
        check(VecGetArrayRead(solution_, &solution_values), "accessing solution");
        out.resize(cells_);
        for (std::size_t i = 0; i < cells_; ++i)
            out[i] =
                static_cast<Real>(PetscRealPart(solution_values[i]));
        check(VecRestoreArrayRead(solution_, &solution_values), "releasing solution");

        return static_cast<Index>(iterations);
    }

    [[nodiscard]] double factorisation_seconds() const noexcept {
        return factorisation_seconds_;
    }

    [[nodiscard]] KSPConvergedReason last_reason() const noexcept {
        return last_reason_;
    }

private:
    static void copy_csr(
        const gss::CsrMatrix& matrix,
        std::vector<PetscInt>& row,
        std::vector<PetscInt>& column,
        std::vector<PetscScalar>& value) {

        const auto r = matrix.row_start();
        const auto c = matrix.column();
        const auto v = matrix.value();

        row.assign(r.begin(), r.end());
        column.assign(c.begin(), c.end());
        value.assign(v.begin(), v.end());
    }

    static void check(PetscErrorCode code, const char* operation) {
        if (code != PETSC_SUCCESS)
            throw std::runtime_error(
                std::string{"PreconditionedNewtonSolver: "} +
                operation + " failed");
    }

    std::size_t cells_{0U};
    std::vector<PetscInt> a_row_, a_col_, j_row_, j_col_;
    std::vector<PetscScalar> a_val_, j_val_, j_base_;
    std::vector<std::size_t> b_to_j_;
    Mat a_matrix_{nullptr};
    Mat j_matrix_{nullptr};
    KSP ksp_{nullptr};
    KSPConvergedReason last_reason_{KSP_CONVERGED_ITERATING};
    Vec solution_{nullptr};
    Vec rhs_{nullptr};
    double factorisation_seconds_{0.0};
};

struct SpectralMeasurement {
    Real rho{std::numeric_limits<Real>::quiet_NaN()};
    std::vector<Real> b_values;
    Index eps_iterations{0};
    Index applications{0};
    Index converged_eigenpairs{0};
    Real max_relative_error{std::numeric_limits<Real>::quiet_NaN()};
    double seconds{0.0};
    double jacobian_seconds{0.0};
};

class SpectralEvaluator final {
public:
    SpectralEvaluator(
        const gss::CsrMatrix& a,
        const gss::CsrMatrix& b_initial,
        const Settings& settings)
        : operator_(
              a,
              b_initial,
              spectral::OperatorMode::Picard,
              Real{1},
              true),
          settings_(settings) {}

    template<class Source, class JacobianAssembler>
    [[nodiscard]] SpectralMeasurement evaluate(
        const ps::CaseWorkspace& workspace,
        std::span<const Real> state,
        const Source& source,
        JacobianAssembler& jacobian_assembler) {

        if (state.size() != workspace.jacobian_pattern.cells)
            throw std::invalid_argument(
                "SpectralEvaluator: state/workspace size mismatch");

        SpectralMeasurement measurement;

        const auto jacobian_start = std::chrono::steady_clock::now();
        measurement.b_values =
            jacobian_assembler.values(state, source);
        measurement.jacobian_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - jacobian_start).count();

        operator_.update_b_values(measurement.b_values);
        operator_.set_mode(spectral::OperatorMode::Picard);

        const auto spectrum =
            spectral::leading_spectrum(
                operator_,
                settings_.requested_eigenpairs,
                settings_.eps_tolerance,
                settings_.eps_max_iterations);

        measurement.rho = spectrum.rho;
        measurement.eps_iterations =
            static_cast<Index>(spectrum.eps_iterations);
        measurement.applications =
            static_cast<Index>(spectrum.applications);
        measurement.converged_eigenpairs =
            static_cast<Index>(spectrum.converged_eigenpairs);
        measurement.max_relative_error =
            static_cast<Real>(spectrum.max_relative_error);
        measurement.seconds = spectrum.eigensolve_seconds;

        return measurement;
    }

    [[nodiscard]] double factorisation_seconds() const noexcept {
        return operator_.factorisation_seconds();
    }

private:
    spectral::MatrixFreePicardOperator operator_;
    const Settings& settings_;
};

template<class Linear>
[[nodiscard]] DenseVector picard_step(
    const ps::CaseWorkspace& workspace,
    Linear& linear,
    std::span<const Real> nonlinear,
    std::span<const Real> state,
    Real omega) {
    DenseVector rhs(nonlinear.size(), Real{0});

    if (rhs.size() != workspace.poisson.boundary_rhs.size())
        throw std::runtime_error("picard_step: RHS size mismatch");

    for (std::size_t i = 0; i < rhs.size(); ++i)
        rhs[i] = nonlinear[i] + workspace.poisson.boundary_rhs[i];

    const auto solved = linear.solve(rhs);

    DenseVector next(state.begin(), state.end());
    for (std::size_t i = 0; i < next.size(); ++i)
        next[i] =
            (Real{1} - omega) * state[i] +
            omega * solved.solution[i];

    return next;
}

template<class Source, class JacobianAssembler>
[[nodiscard]] std::pair<DenseVector, Index> newton_step(
    PreconditionedNewtonSolver& newton,
    const Source& source,
    std::span<const Real> state,
    Real damping,
    JacobianAssembler& jacobian_assembler,
    std::span<const Real> r,
    const std::vector<Real>* preassembled_b_values = nullptr, double* assembly_seconds=nullptr) {

    std::vector<Real> assembled_b;
    if (preassembled_b_values == nullptr) {
        const auto begin=std::chrono::steady_clock::now();
        assembled_b = jacobian_assembler.values(state, source);
        if(assembly_seconds) *assembly_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        preassembled_b_values = &assembled_b;
    }
    newton.update_source_jacobian_values(*preassembled_b_values);

    DenseVector rhs(r.size(), Real{0});
    for (std::size_t i = 0; i < r.size(); ++i)
        rhs[i] = -r[i];

    DenseVector correction;
    const Index krylov_iterations =
        newton.solve(rhs, correction);

    DenseVector next(state.begin(), state.end());
    for (std::size_t i = 0; i < next.size(); ++i)
        next[i] += damping * correction[i];

    return {std::move(next), krylov_iterations};
}


struct RandomCase {
    std::size_t case_id{0U};
    Real alpha{Real{0}};
    Real theta_degrees{Real{0}};
    Real rho_target{Real{0}};
    Real rho_calibration_c1{std::numeric_limits<Real>::quiet_NaN()};
    Real c{std::numeric_limits<Real>::quiet_NaN()};
    std::string calibration_error;
};

[[nodiscard]] inline std::vector<RandomCase>
generate_random_cases(const Settings& settings) {
    std::mt19937_64 generator(settings.sampling_seed);
    std::uniform_real_distribution<double> alpha_dist(
        static_cast<double>(settings.alpha_min),
        static_cast<double>(settings.alpha_max));
    std::uniform_real_distribution<double> theta_dist(
        static_cast<double>(settings.theta_min_degrees),
        static_cast<double>(settings.theta_max_degrees));
    std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

    std::vector<RandomCase> cases;
    cases.reserve(settings.sample_count);

    // The stratum is assigned deterministically by case_id modulo five.
    // Samples are distributed evenly, with any remainder in the first strata.
    // Alpha, theta and rho_target are generated before any nonlinear solve,
    // using the same generator and the same call order in all four executables.
    constexpr std::size_t strata = 5U;
    for (std::size_t q = 0; q < settings.sample_count; ++q) {
        const std::size_t bin = q % strata;
        const Real alpha = static_cast<Real>(alpha_dist(generator));
        const Real theta_degrees = static_cast<Real>(theta_dist(generator));
        const double u = unit_dist(generator);
        const double lo = static_cast<double>(settings.rho_strata_edges[bin]);
        const double hi = static_cast<double>(settings.rho_strata_edges[bin + 1U]);
        const Real rho_target = static_cast<Real>(lo + (hi - lo) * u);

        cases.push_back(RandomCase{
            q + 1U,
            alpha,
            theta_degrees,
            rho_target,
            std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN(),
            {}});
    }
    return cases;
}

// Only Fadeev data change between samples; geometry and A remain invariant.
inline void refresh_case_data(ps::CaseWorkspace& workspace,
                              const gss::Fadeev& fadeev) {
    const auto& topology = workspace.mesh.topology();
    for (std::size_t i = 0; i < workspace.exact.size(); ++i)
        workspace.exact[i] = fadeev.psi(topology.coord_x[i], topology.coord_y[i]);

    const gss::BoundaryCondition condition = gss::Dirichlet{
        [&fadeev](const gss::BoundaryPoint2D& point) {
            return fadeev.psi(point.x, point.y);
        }};
    const auto boundary = gss::BoundaryApplication::evaluate(
        workspace.mesh.geometry().boundary_faces(), condition);
    auto& rhs = workspace.poisson.boundary_rhs;
    std::fill(rhs.begin(), rhs.end(), Real{0});
    // Same face order and arithmetic as PoissonAssembler::assemble, with
    // fixed Dirichlet coefficients alpha=1, beta=0 and diffusivity=1.
    for (const auto& face : boundary) {
        const auto owner = static_cast<std::size_t>(face.geometry.owner_cell);
        const Real d = face.coefficients.alpha * face.geometry.perpendicular_distance
                     + face.coefficients.beta;
        if (owner >= rhs.size() || !(d > Real{0}))
            throw std::runtime_error("invalid invariant Dirichlet boundary geometry");
        rhs[owner] += Real{1} * face.geometry.measure * face.coefficients.gamma / d;
    }
    workspace.set_nonlinearity_scale(fadeev.nonlinearity_scale());
}

template<class GradientMethod>
void calibrate_case_scales(
    const Settings& settings,
    std::vector<RandomCase>& random_cases,
    const gss::Config& base_config,
    ps::StudySettings study_settings,
    const GradientMethod& method) {

    const Real delta = base_config.solver.source_regularisation_delta;
    const auto mesh_config =
        ps::make_mesh_config(
            base_config.mesh,
            "cartesian",
            settings.mesh_seed,
            settings.nominal_nx,
            settings.nominal_ny);

    const gss::Fadeev calibration_fadeev{
        settings.calibration_alpha,
        settings.calibration_theta_degrees,
        Real{1}};

    auto seed_workspace =
        ps::make_workspace(
            mesh_config,
            calibration_fadeev,
            study_settings.jacobian_stencil_radius);

    const DenseVector seed_reference = seed_workspace.exact;
    const auto seed_base_source =
        evaluate_base_source(
            seed_workspace,
            seed_reference,
            calibration_fadeev,
            method,
            delta,
            study_settings);
    const auto seed_compensator =
        make_compensator(
            seed_workspace,
            seed_reference,
            seed_base_source,
            Real{1});
    const ModifiedSource<GradientMethod> seed_source{
        seed_workspace,
        calibration_fadeev,
        method,
        delta,
        study_settings,
        Real{1},
        seed_compensator};

    SourceJacobianAssembler<GradientMethod> jacobian_assembler{
        seed_workspace,
        method,
        delta,
        study_settings,
        settings.source_derivative};

    const auto seed_b_values =
        jacobian_assembler.values(seed_reference, seed_source);
    const auto seed_b =
        sparse::to_csr(seed_workspace.jacobian_pattern, seed_b_values);

    SpectralEvaluator evaluator{
        seed_workspace.poisson.matrix,
        seed_b,
        settings};

    std::cout
        << "Calibrating c from rho_target on canonical Cartesian/WLSQR "
        << "reference states...\n";

    for (auto& sample : random_cases) {
        try {
            ps::StudySettings local_study = study_settings;
            local_study.alpha = sample.alpha;
            local_study.theta_degrees = sample.theta_degrees;

            const gss::Fadeev fadeev{
                sample.alpha,
                sample.theta_degrees,
                Real{1}};

            auto& workspace = seed_workspace;
            refresh_case_data(workspace, fadeev);

            const DenseVector reference_state = workspace.exact;
            const auto base_source_at_reference =
                evaluate_base_source(
                    workspace,
                    reference_state,
                    fadeev,
                    method,
                    delta,
                    local_study);
            const auto compensator_c1 =
                make_compensator(
                    workspace,
                    reference_state,
                    base_source_at_reference,
                    Real{1});
            const ModifiedSource<GradientMethod> source_c1{
                workspace,
                fadeev,
                method,
                delta,
                local_study,
                Real{1},
                compensator_c1};

            const auto spectrum_c1 =
                evaluator.evaluate(
                    workspace,
                    reference_state,
                    source_c1,
                    jacobian_assembler);

            sample.rho_calibration_c1 = spectrum_c1.rho;

            if (!std::isfinite(static_cast<double>(sample.rho_calibration_c1)) ||
                !(sample.rho_calibration_c1 > Real{0}))
                throw std::runtime_error(
                    "non-positive or non-finite rho(G) at c=1");

            sample.c = sample.rho_target / sample.rho_calibration_c1;
            if (!std::isfinite(static_cast<double>(sample.c)) ||
                !(sample.c > Real{0}))
                throw std::runtime_error(
                    "calibrated c is non-positive or non-finite");

            sample.calibration_error.clear();
        } catch (const std::exception& error) {
            sample.rho_calibration_c1 =
                std::numeric_limits<Real>::quiet_NaN();
            sample.c = std::numeric_limits<Real>::quiet_NaN();
            sample.calibration_error = error.what();
        }

        if (settings.progress_every > 0U &&
            sample.case_id % settings.progress_every == 0U) {
            std::cout
                << "calibration case " << sample.case_id
                << "/" << settings.sample_count
                << " rho_target=" << sample.rho_target
                << " rho_c1=" << sample.rho_calibration_c1
                << " c=" << sample.c << '\n';
        }
    }
}

// Semicolon-separated entries keep the field valid in the comma-separated CSV.
[[nodiscard]] inline std::string linear_reason_text(
    const std::map<int, std::size_t>& counts) {
    std::ostringstream out;
    bool first = true;
    for (const auto& [reason, count] : counts) {
        if (!first) out << ';';
        first = false;
        const char* name = KSPConvergedReasons[reason];
        out << reason << ':' << (name ? name : "unknown") << ':' << count;
    }
    return out.str();
}

struct RunSummary {
    bool converged{false};
    DenseVector final_state;
    std::string termination{"not_run"};

    std::size_t outer_iterations{0U};
    std::size_t picard_steps{0U};
    std::size_t newton_steps{0U};
    std::size_t rho_evaluations{0U};
    std::size_t solver_switches{0U};

    std::size_t picard_linear_solve_calls{0U};
    std::size_t newton_linear_solve_calls{0U};
    std::size_t spectral_A_solve_calls{0U};
    std::size_t newton_linear_hits_50{0U};
    std::size_t linear_iteration_limits{0U};
    std::size_t linear_residual_divergences{0U};
    std::size_t linear_other_failures{0U};
    std::map<int, std::size_t> linear_reason_counts;
    Index krylov_iterations_total{0};

    Index eps_iterations_total{0};
    Index eps_applications_total{0};

    Real initial_residual{std::numeric_limits<Real>::quiet_NaN()};
    Real final_residual{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_target{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_reference{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_initial{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_last_decision{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_final_diagnostic{std::numeric_limits<Real>::quiet_NaN()};

    Real initial_reference_error{std::numeric_limits<Real>::quiet_NaN()};
    Real final_reference_error{std::numeric_limits<Real>::quiet_NaN()};

    bool rho_predicts_picard_convergence{false};
    bool prediction_failure_residual{false};
    bool prediction_failure_reference{false};
    bool converged_to_reference{false};
    bool unexpected_picard_reference_convergence{false};

    double seconds_algorithmic{0.0};
    double seconds_spectral_strategy{0.0};
    double seconds_jacobian_assembly{0.0};
    std::size_t regularisation_floor_active_iterations{0};
    double seconds_picard{0.0};
    double seconds_newton{0.0};
};

struct CaseMetric {
    std::string method;
    std::string mesh;
    Real rho_target{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_reference{std::numeric_limits<Real>::quiet_NaN()};
    Real rho_initial{std::numeric_limits<Real>::quiet_NaN()};
    bool converged{false};
    bool case_exception{false};
    bool converged_to_reference{false};
    bool prediction_failure_residual{false};
    bool prediction_failure_reference{false};
    bool unexpected_picard_reference_convergence{false};
    std::size_t linear_solve_calls{0U};
    std::size_t linear_hits_50{0U};
    std::size_t rho_evaluations{0U};
    std::size_t solver_switches{0U};
    Index krylov_iterations{0};
    double seconds_algorithmic{0.0};
    double seconds_case_total{0.0};
    double seconds_initial_rho_diagnostic{0.0};
    double seconds_final_rho_diagnostic{0.0};
};

[[nodiscard]] inline double quantile(
    std::vector<double> values,
    double probability) {

    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();

    std::sort(values.begin(), values.end());
    const double position =
        probability * static_cast<double>(values.size() - 1U);
    const auto lo = static_cast<std::size_t>(std::floor(position));
    const auto hi = static_cast<std::size_t>(std::ceil(position));

    if (lo == hi)
        return values[lo];

    const double w = position - static_cast<double>(lo);
    return (1.0 - w) * values[lo] + w * values[hi];
}

[[nodiscard]] inline double mean(const std::vector<double>& values) {
    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();
    return std::accumulate(values.begin(), values.end(), 0.0) /
           static_cast<double>(values.size());
}

inline void write_summary_row(
    std::ofstream& out,
    std::string_view scope,
    std::string_view method,
    std::string_view mesh,
    const std::vector<CaseMetric>& records) {

    if (records.empty())
        return;

    std::size_t converged = 0U;
    std::size_t converged_ref = 0U;
    std::size_t failed = 0U;
    std::size_t case_exceptions = 0U;
    std::size_t pred_fail_res = 0U;
    std::size_t pred_fail_ref = 0U;
    std::size_t unexpected = 0U;
    std::size_t linear_calls = 0U;
    std::size_t linear_hits_50 = 0U;
    std::size_t rho_evals = 0U;
    std::size_t switches = 0U;
    long long krylov = 0;

    std::vector<double> algorithmic_times;
    std::vector<double> total_times;
    algorithmic_times.reserve(records.size());
    total_times.reserve(records.size());

    for (const auto& r : records) {
        converged += r.converged ? 1U : 0U;
        converged_ref += r.converged_to_reference ? 1U : 0U;
        failed += r.converged ? 0U : 1U;
        case_exceptions += r.case_exception ? 1U : 0U;
        pred_fail_res += r.prediction_failure_residual ? 1U : 0U;
        pred_fail_ref += r.prediction_failure_reference ? 1U : 0U;
        unexpected += r.unexpected_picard_reference_convergence ? 1U : 0U;
        linear_calls += r.linear_solve_calls;
        linear_hits_50 += r.linear_hits_50;
        rho_evals += r.rho_evaluations;
        switches += r.solver_switches;
        krylov += static_cast<long long>(r.krylov_iterations);
        algorithmic_times.push_back(r.seconds_algorithmic);
        total_times.push_back(r.seconds_case_total);
    }

    out
        << scope << ','
        << method << ','
        << mesh << ','
        << records.size() << ','
        << converged << ','
        << failed << ','
        << case_exceptions << ','
        << converged_ref << ','
        << pred_fail_res << ','
        << pred_fail_ref << ','
        << unexpected << ','
        << linear_calls << ','
        << static_cast<double>(linear_calls) / static_cast<double>(records.size()) << ','
        << linear_hits_50 << ','
        << rho_evals << ','
        << switches << ','
        << krylov << ','
        << mean(algorithmic_times) << ','
        << quantile(algorithmic_times, 0.50) << ','
        << quantile(algorithmic_times, 0.90) << ','
        << quantile(algorithmic_times, 0.95) << ','
        << *std::max_element(algorithmic_times.begin(), algorithmic_times.end()) << ','
        << mean(total_times) << ','
        << quantile(total_times, 0.50) << ','
        << quantile(total_times, 0.95)
        << '\n';
}

inline void write_summary(
    const std::filesystem::path& result_dir,
    const std::vector<CaseMetric>& all_records, const Settings& settings) {

    std::ofstream out(result_dir / "summary.csv");
    if (!out)
        throw std::runtime_error("could not create summary.csv");

    out << std::scientific << std::setprecision(17);
    out
        << "scope,method,mesh,cases,converged,not_converged,case_exceptions,"
           "converged_to_reference,prediction_failure_residual,"
           "prediction_failure_reference,unexpected_picard_reference_convergence,"
           "linear_solve_calls,mean_linear_solve_calls_per_case,"
           "linear_hits_50,rho_evaluations,solver_switches,"
           "krylov_iterations_total,mean_algorithmic_seconds,"
           "median_algorithmic_seconds,p90_algorithmic_seconds,"
           "p95_algorithmic_seconds,max_algorithmic_seconds,"
           "mean_case_total_seconds,median_case_total_seconds,"
           "p95_case_total_seconds\n";

    write_summary_row(out, "overall", "all", "all", all_records);

    const auto& methods=settings.methods;
    const auto& meshes=settings.mesh_families;

    for (const auto& method : methods) {
        for (const auto& mesh : meshes) {
            std::vector<CaseMetric> selected;
            for (const auto& r : all_records)
                if (r.method == method && r.mesh == mesh)
                    selected.push_back(r);
            write_summary_row(out, "method_mesh", method, mesh, selected);
        }
    }
}

inline void write_rho_bins(
    const std::filesystem::path& result_dir,
    const std::vector<CaseMetric>& records, bool by_target=false) {

    const std::array<double, 8> edges{
        0.0, 0.70, 0.90, 0.98, 1.02, 1.10, 1.30,
        std::numeric_limits<double>::infinity()};

    std::ofstream out(result_dir / "rho_bins.csv");
    if (!out)
        throw std::runtime_error("could not create rho_bins.csv");

    out << std::scientific << std::setprecision(17);
    out
        << "rho_min,rho_max,cases,converged,converged_to_reference,"
           "prediction_failure_residual,prediction_failure_reference,"
           "mean_linear_solve_calls,mean_algorithmic_seconds,"
           "median_algorithmic_seconds,p95_algorithmic_seconds\n";

    for (std::size_t b = 0; b + 1U < edges.size(); ++b) {
        const double lo = edges[b];
        const double hi = edges[b + 1U];

        std::vector<const CaseMetric*> selected;
        for (const auto& r : records) {
            const double rho = static_cast<double>(by_target ? r.rho_target : r.rho_reference);
            if (!std::isfinite(rho))
                continue;
            if (rho >= lo && rho < hi)
                selected.push_back(&r);
        }

        if (selected.empty())
            continue;

        std::size_t conv = 0U;
        std::size_t conv_ref = 0U;
        std::size_t pfr = 0U;
        std::size_t pff = 0U;
        std::size_t linear_calls = 0U;
        std::vector<double> times;

        for (const auto* r : selected) {
            conv += r->converged ? 1U : 0U;
            conv_ref += r->converged_to_reference ? 1U : 0U;
            pfr += r->prediction_failure_residual ? 1U : 0U;
            pff += r->prediction_failure_reference ? 1U : 0U;
            linear_calls += r->linear_solve_calls;
            times.push_back(r->seconds_algorithmic);
        }

        out
            << lo << ',' << hi << ','
            << selected.size() << ','
            << conv << ','
            << conv_ref << ','
            << pfr << ','
            << pff << ','
            << static_cast<double>(linear_calls) /
               static_cast<double>(selected.size()) << ','
            << mean(times) << ','
            << quantile(times, 0.50) << ','
            << quantile(times, 0.95)
            << '\n';
    }
}

template<class Linear, class Source, class JacobianAssembler>
[[nodiscard]] RunSummary solve_case(
    Strategy strategy,
    const Settings& settings,
    const ps::CaseWorkspace& workspace,
    Linear& picard_linear,
    PreconditionedNewtonSolver* newton,
    SpectralEvaluator& spectral_evaluator,
    const Source& source,
    JacobianAssembler& jacobian_assembler,
    DenseVector state,
    std::span<const Real> reference_state,
    const SpectralMeasurement& reference_spectral,
    const SpectralMeasurement& initial_spectral,
    double initial_spectral_seconds,
    std::ofstream& switches,
    std::size_t case_id,
    std::string_view method_name,
    std::string_view mesh_name,
    Real alpha,
    Real theta_degrees,
    Real c,
    Real rho_target) {

    RunSummary summary;
    summary.final_state = state;
    summary.rho_target = rho_target;
    summary.rho_reference = reference_spectral.rho;
    summary.rho_initial = initial_spectral.rho;
    summary.rho_predicts_picard_convergence =
        std::isfinite(static_cast<double>(reference_spectral.rho)) &&
        reference_spectral.rho < settings.picard_if_rho_below;
    summary.initial_reference_error =
        ps::relative_difference(state, reference_state);

    const auto solve_start = std::chrono::steady_clock::now();
    auto nonlinear = source(state);
    auto r = residual_from_source(workspace, nonlinear, state);
    const Real initial_norm = ps::norm2(r);
    summary.initial_residual = initial_norm;
    summary.final_residual = initial_norm;

    const Real threshold =
        settings.tolerance_absolute +
        settings.tolerance_relative * std::max(initial_norm, Real{1});

    const Real divergence_limit =
        settings.divergence_factor *
        std::max(initial_norm, Real{1});

    bool fixed_picard = strategy == Strategy::PicardOnly;

    std::vector<Real> reusable_b_values;
    bool have_reusable_b_values = false;
    double reusable_b_seconds=initial_spectral.jacobian_seconds;

    if (strategy == Strategy::MixedRhoInitial ||
        strategy == Strategy::MixedRhoEveryIteration) {

        summary.rho_last_decision = initial_spectral.rho;
        ++summary.rho_evaluations;
        summary.seconds_spectral_strategy += initial_spectral_seconds;
        summary.spectral_A_solve_calls +=
            static_cast<std::size_t>(initial_spectral.applications);
        summary.eps_iterations_total += initial_spectral.eps_iterations;
        summary.eps_applications_total += initial_spectral.applications;

        fixed_picard =
            initial_spectral.rho < settings.picard_if_rho_below;

        reusable_b_values = initial_spectral.b_values;
        have_reusable_b_values = !fixed_picard;
    }

    bool previous_method_valid = false;
    bool previous_was_picard = false;

    for (std::size_t k = 0; k < settings.max_iterations; ++k) {
        if (k != 0U) {
            nonlinear = source(state);
            r = residual_from_source(workspace, nonlinear, state);
        }
        const Real rn = ps::norm2(r);
        summary.final_residual = rn;

        if (!std::isfinite(static_cast<double>(rn))) {
            summary.termination = "nonfinite_residual";
            break;
        }
        if (rn <= threshold) {
            summary.converged = true;
            summary.termination = "converged";
            break;
        }
        if (rn > divergence_limit) {
            summary.termination = "diverged";
            break;
        }

        if(settings.floor_activity && source.floor_active()) ++summary.regularisation_floor_active_iterations;
        bool use_picard = fixed_picard;
        Real rho_used = std::numeric_limits<Real>::quiet_NaN();
        const std::vector<Real>* b_for_newton = nullptr;

        if (strategy == Strategy::NewtonOnly) {
            use_picard = false;
        } else if (strategy == Strategy::MixedRhoEveryIteration) {
            if (k == 0U) {
                rho_used = initial_spectral.rho;
            } else {
                const auto t = std::chrono::steady_clock::now();
                const auto spectral =
                    spectral_evaluator.evaluate(
                        workspace, state, source, jacobian_assembler);
                const double wall =
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t).count();

                rho_used = spectral.rho;
                summary.rho_last_decision = rho_used;
                ++summary.rho_evaluations;
                summary.seconds_spectral_strategy += wall;
                summary.spectral_A_solve_calls +=
                    static_cast<std::size_t>(spectral.applications);
                summary.eps_iterations_total += spectral.eps_iterations;
                summary.eps_applications_total += spectral.applications;

                reusable_b_values = spectral.b_values;
                reusable_b_seconds = spectral.jacobian_seconds;
                have_reusable_b_values =
                    !(rho_used < settings.picard_if_rho_below);
            }

            use_picard = rho_used < settings.picard_if_rho_below;
        } else if (strategy == Strategy::MixedRhoInitial) {
            rho_used = initial_spectral.rho;
        }

        if (previous_method_valid &&
            previous_was_picard != use_picard) {

            ++summary.solver_switches;
            switches
                << case_id << ','
                << method_name << ','
                << mesh_name << ','
                << alpha << ','
                << theta_degrees << ','
                << rho_target << ','
                << c << ','
                << k << ','
                << rn << ','
                << rho_used << ','
                << (previous_was_picard ? "picard" : "newton") << ','
                << (use_picard ? "picard" : "newton")
                << '\n';
            switches.flush();
        }

        previous_method_valid = true;
        previous_was_picard = use_picard;

        try {
            if (use_picard) {
                const auto t = std::chrono::steady_clock::now();
                state =
                    picard_step(
                        workspace,
                        picard_linear,
                        nonlinear,
                        state,
                        settings.picard_omega);
                summary.seconds_picard +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t).count();

                ++summary.picard_steps;
                ++summary.picard_linear_solve_calls;
                have_reusable_b_values = false;
            } else {
                if (newton == nullptr)
                    throw std::runtime_error(
                        "Newton requested but Newton solver was not constructed");

                if (have_reusable_b_values)
                    b_for_newton = &reusable_b_values;
                else if (k == 0U && strategy == Strategy::NewtonOnly && !initial_spectral.b_values.empty())
                    b_for_newton = &initial_spectral.b_values;

                if(settings.jacobian_timing && b_for_newton)
                    summary.seconds_jacobian_assembly += have_reusable_b_values ? reusable_b_seconds : initial_spectral.jacobian_seconds;
                const auto t = std::chrono::steady_clock::now();
                auto [next, krylov] =
                    newton_step(
                        *newton,
                        source,
                        state,
                        settings.newton_damping,
                        jacobian_assembler,
                        r,
                        b_for_newton, settings.jacobian_timing ? &summary.seconds_jacobian_assembly : nullptr);
                summary.seconds_newton +=
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - t).count();

                // Observation only: keep using the correction, including when
                // the deliberately small linear iteration budget is exhausted.
                const auto reason = newton->last_reason();
                ++summary.linear_reason_counts[static_cast<int>(reason)];
                if (reason == KSP_DIVERGED_ITS)
                    ++summary.linear_iteration_limits;
                else if (reason == KSP_DIVERGED_DTOL)
                    ++summary.linear_residual_divergences;
                else if (reason < 0)
                    ++summary.linear_other_failures;

                state = std::move(next);
                summary.krylov_iterations_total += krylov;
                ++summary.newton_steps;
                ++summary.newton_linear_solve_calls;
                if (krylov >= settings.krylov_max_iterations)
                    ++summary.newton_linear_hits_50;
                have_reusable_b_values = false;
            }
        } catch (const std::exception& error) {
            summary.termination =
                std::string{"step_failure:"} + error.what();
            break;
        }

        summary.outer_iterations = k + 1U;
    }

    if (!summary.converged &&
        summary.termination == "not_run")
        summary.termination = "max_iterations";

    summary.seconds_algorithmic =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() - solve_start).count();

    // For mixed strategies the initial spectral computation is part of the
    // strategy itself, although it was evaluated immediately before entry.
    if (strategy == Strategy::MixedRhoInitial ||
        strategy == Strategy::MixedRhoEveryIteration)
        summary.seconds_algorithmic += initial_spectral_seconds;

    // Charge the reused initial B to Newton, not to the diagnostic. It is
    // required by the algorithm even though the diagnostic assembled it first.
    if (strategy == Strategy::NewtonOnly && summary.newton_steps > 0) {
        summary.seconds_algorithmic += initial_spectral.jacobian_seconds;
        summary.seconds_newton += initial_spectral.jacobian_seconds;
    }

    summary.final_state = state;
    summary.final_reference_error =
        ps::relative_difference(state, reference_state);

    summary.converged_to_reference =
        summary.converged &&
        std::isfinite(static_cast<double>(summary.final_reference_error)) &&
        summary.final_reference_error <= settings.reference_relative_tolerance;

    summary.prediction_failure_residual =
        summary.rho_predicts_picard_convergence &&
        !summary.converged;

    summary.prediction_failure_reference =
        summary.rho_predicts_picard_convergence &&
        !summary.converged_to_reference;

    summary.unexpected_picard_reference_convergence =
        strategy == Strategy::PicardOnly &&
        !summary.rho_predicts_picard_convergence &&
        summary.converged_to_reference;

    return summary;
}


inline std::string csv_safe(std::string text) {
    for (char& ch : text) {
        if (ch == ',' || ch == '\n' || ch == '\r')
            ch = ';';
    }
    return text;
}

template<class Sample>
void write_failed_case(
    std::ofstream& cases,
    const Sample& sample,
    std::string_view method_name,
    std::string_view mesh_name,
    Index actual_cells,
    Strategy strategy,
    std::string_view initial_guess_kind,
    const std::string& reason,
    double seconds_case_total) {

    const Real nan = std::numeric_limits<Real>::quiet_NaN();

    cases
        << sample.case_id << ','
        << sample.alpha << ','
        << sample.theta_degrees << ','
        << sample.rho_target << ','
        << sample.rho_calibration_c1 << ','
        << sample.c << ','
        << method_name << ','
        << mesh_name << ','
        << actual_cells << ','
        << strategy_name(strategy) << ','
        << initial_guess_kind << ','
        << nan << ','                             // rho_reference
        << nan << ','                             // rho_initial
        << nan << ','                             // rho_last_decision
        << nan << ','                             // rho_final_diagnostic
        << 0 << ','                               // rho predicts Picard
        << 0 << ','                               // converged
        << 0 << ','                               // converged_to_reference
        << 0 << ','                               // prediction_failure_residual
        << 0 << ','                               // prediction_failure_reference
        << 0 << ','                               // unexpected convergence
        << "case_exception:" << csv_safe(reason) << ','
        << 0 << ','                               // outer_iterations
        << 0 << ','                               // picard_steps
        << 0 << ','                               // newton_steps
        << 0 << ','                               // picard linear calls
        << 0 << ','                               // newton linear calls
        << 0 << ','                               // strategy spectral A calls
        << 0 << ','                               // diagnostic initial A calls
        << 0 << ','                               // diagnostic reference A calls
        << 0 << ','                               // diagnostic final A calls
        << 0 << ','                               // strategy linear calls
        << 0 << ','                               // total observed linear calls
        << 0 << ','                               // Krylov iterations
        << 0 << ','                               // Newton hits 50
        << 0 << ','                               // rho evaluations
        << 0 << ','                               // solver switches
        << 0 << ','                               // EPS iterations
        << 0 << ','                               // EPS applications
        << nan << ','                             // initial residual
        << nan << ','                             // final residual
        << nan << ','                             // initial reference error
        << nan << ','                             // final reference error
        << 0.0 << ','                             // algorithmic time
        << 0.0 << ','                             // reference rho diagnostic
        << 0.0 << ','                             // initial rho diagnostic
        << 0.0 << ','                             // spectral strategy
        << 0.0 << ','                             // Picard
        << 0.0 << ','                             // Newton
        << 0.0 << ','                             // final rho diagnostic
        << 0.0 << ','                             // case setup
        << seconds_case_total
        << ",0,0,0,,nan,nan" // case failed before a valid accumulated observation could be returned
        << '\n';

    cases.flush();
}

template<class GradientMethod>
void run_method(
    Strategy strategy,
    std::string_view method_name,
    const GradientMethod& method,
    const Settings& settings,
    const std::vector<RandomCase>& random_cases,
    const gss::Config& base_config,
    ps::StudySettings study_settings,
    std::ofstream& cases,
    std::ofstream& switches,
    std::ofstream& setup,
    std::vector<CaseMetric>& metrics,
    std::size_t& algorithmic_nonconvergences,
    std::size_t& prediction_failures_residual,
    std::size_t& prediction_failures_reference,
    std::size_t& wrong_equilibria,
    std::size_t& case_failures,
    std::size_t& setup_failures) {

    const Real delta =
        base_config.solver.source_regularisation_delta;

    for (const auto& mesh_name : settings.mesh_families) {
        try {
            const auto mesh_config =
                ps::make_mesh_config(
                    base_config.mesh,
                    mesh_name,
                    settings.mesh_seed,
                    settings.nominal_nx,
                    settings.nominal_ny);

            // One invariant A-factorisation per mesh/method campaign.  The
            // Fadeev parameters here only provide the boundary data needed to
            // construct a workspace; A and its adjacency pattern are mesh-only.
            const gss::Fadeev calibration_fadeev{
                settings.calibration_alpha,
                settings.calibration_theta_degrees,
                Real{1}};

            auto solver_workspace =
                ps::make_workspace(
                    mesh_config,
                    calibration_fadeev,
                    study_settings.jacobian_stencil_radius);

            const auto actual_cells =
                static_cast<Index>(solver_workspace.exact.size());

            auto picard_linear =
                ps::make_picard_controller(solver_workspace);

            const DenseVector calibration_reference =
                solver_workspace.exact;

            const auto calibration_base_source =
                evaluate_base_source(
                    solver_workspace,
                    calibration_reference,
                    calibration_fadeev,
                    method,
                    delta,
                    study_settings);

            const auto calibration_compensator =
                make_compensator(
                    solver_workspace,
                    calibration_reference,
                    calibration_base_source,
                    Real{1});

            const ModifiedSource<GradientMethod> calibration_source{
                solver_workspace,
                calibration_fadeev,
                method,
                delta,
                study_settings,
                Real{1},
                calibration_compensator};

            SourceJacobianAssembler<GradientMethod> jacobian_assembler{
                solver_workspace,
                method,
                delta,
                study_settings,
                settings.source_derivative};

            const auto calibration_b_values =
                jacobian_assembler.values(
                    calibration_reference,
                    calibration_source);

            const auto calibration_b =
                sparse::to_csr(
                    solver_workspace.jacobian_pattern,
                    calibration_b_values);

            SpectralEvaluator spectral_evaluator{
                solver_workspace.poisson.matrix,
                calibration_b,
                settings};

            const auto calibration_j =
                sparse::subtract_from(
                    solver_workspace.poisson.matrix,
                    solver_workspace.jacobian_pattern,
                    calibration_b_values);

            if (solver_workspace.jacobian_pattern.radius == 1U)
                require_matrix_ring_pattern(
                    calibration_j,
                    solver_workspace.mesh.topology(),
                    1U,
                    "Newton J");

            std::unique_ptr<PreconditionedNewtonSolver> newton;
            if (strategy != Strategy::PicardOnly) {
                newton =
                    std::make_unique<PreconditionedNewtonSolver>(
                        solver_workspace.poisson.matrix,
                        calibration_j,
                        solver_workspace.jacobian_pattern,
                        true,
                        settings.krylov_relative_tolerance,
                        settings.krylov_max_iterations);
            }

            const auto calibration_spectrum =
                spectral_evaluator.evaluate(
                    solver_workspace,
                    calibration_reference,
                    calibration_source,
                    jacobian_assembler);

            setup
                << method_name << ','
                << mesh_name << ','
                << settings.target_cells << ','
                << actual_cells << ','
                << settings.nominal_nx << ','
                << settings.nominal_ny << ','
                << solver_workspace.jacobian_pattern.colours << ','
                << solver_workspace.jacobian_pattern.nnz() << ','
                << jacobian_assembler.derivative_method() << ','
                << jacobian_assembler.gradient_operator_relative_defect() << ','
                << calibration_spectrum.rho << ','
                << spectral_evaluator.factorisation_seconds() << ','
                << (newton ? newton->factorisation_seconds() : 0.0)
                << '\n';
            setup.flush();

            for (const auto& sample : random_cases) {
                const auto case_wall_start =
                    std::chrono::steady_clock::now();

                const std::string failure_initial_guess_kind =
                    "reference_plus_deterministic_perturbation";

                try {

                if (!sample.calibration_error.empty() ||
                    !std::isfinite(static_cast<double>(sample.c)) ||
                    !(sample.c > Real{0}))
                    throw std::runtime_error(
                        std::string{"rho-target calibration failed: "} +
                        sample.calibration_error);

                ps::StudySettings local_study = study_settings;
                local_study.alpha = sample.alpha;
                local_study.theta_degrees = sample.theta_degrees;

                const gss::Fadeev fadeev{
                    sample.alpha,
                    sample.theta_degrees,
                    Real{1}};

                const auto case_setup_start =
                    std::chrono::steady_clock::now();

                auto& workspace = solver_workspace;
                refresh_case_data(workspace, fadeev);

                if (workspace.exact.size() !=
                    static_cast<std::size_t>(actual_cells))
                    throw std::runtime_error(
                        "case mesh size changed within a fixed mesh campaign");

                const DenseVector reference_state =
                    workspace.exact;

                const auto base_source_at_reference =
                    evaluate_base_source(
                        workspace,
                        reference_state,
                        fadeev,
                        method,
                        delta,
                        local_study);

                const auto compensator =
                    make_compensator(
                        workspace,
                        reference_state,
                        base_source_at_reference,
                        sample.c);

                const ModifiedSource<GradientMethod> source{
                    workspace,
                    fadeev,
                    method,
                    delta,
                    local_study,
                    sample.c,
                    compensator, settings.floor_activity};

                const double case_setup_seconds =
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - case_setup_start).count();

                DenseVector initial =
                    perturbed_state(
                        workspace,
                        reference_state,
                        settings.perturbation_amplitude);

                const char* initial_guess_kind =
                    "reference_plus_deterministic_perturbation";

                const Real initial_reference_error =
                    ps::relative_difference(initial, reference_state);

                SpectralMeasurement reference_spectral, initial_spectral;
                double reference_rho_seconds=0.0, initial_rho_seconds=0.0;
                if(settings.reference_rho) {
                    const auto start=std::chrono::steady_clock::now();
                    reference_spectral=spectral_evaluator.evaluate(workspace,reference_state,source,jacobian_assembler);
                    reference_rho_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                }
                const bool initial_required=strategy==Strategy::MixedRhoInitial || strategy==Strategy::MixedRhoEveryIteration;
                if(settings.initial_rho || initial_required) {
                    const auto start=std::chrono::steady_clock::now();
                    initial_spectral=spectral_evaluator.evaluate(workspace,initial,source,jacobian_assembler);
                    initial_rho_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
                }

                auto run =
                    solve_case(
                        strategy,
                        settings,
                        workspace,
                        picard_linear,
                        newton.get(),
                        spectral_evaluator,
                        source,
                        jacobian_assembler,
                        std::move(initial),
                        reference_state,
                        reference_spectral,
                        initial_spectral,
                        initial_rho_seconds,
                        switches,
                        sample.case_id,
                        method_name,
                        mesh_name,
                        sample.alpha,
                        sample.theta_degrees,
                        sample.c,
                        sample.rho_target);

                const double initial_rho_diagnostic_seconds =
                    settings.initial_rho ? initial_rho_seconds -
                    ((strategy == Strategy::NewtonOnly && run.newton_steps > 0)
                        ? initial_spectral.jacobian_seconds : 0.0) : 0.0;
                if(!settings.initial_rho) run.rho_initial=std::numeric_limits<Real>::quiet_NaN();
                run.initial_reference_error =
                    initial_reference_error;

                double final_rho_seconds = 0.0;
                std::size_t diagnostic_spectral_A_solve_calls = 0U;

                if (strategy == Strategy::PicardOnly &&
                    settings.compute_unexpected_picard_final_rho &&
                    run.converged &&
                    std::isfinite(static_cast<double>(run.rho_reference)) &&
                    run.rho_reference >= settings.picard_if_rho_below) {

                    const auto final_rho_start =
                        std::chrono::steady_clock::now();
                    const auto final_spectral =
                        spectral_evaluator.evaluate(
                            workspace,
                            run.final_state,
                            source,
                            jacobian_assembler);
                    final_rho_seconds =
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - final_rho_start).count();
                    run.rho_final_diagnostic = final_spectral.rho;
                    diagnostic_spectral_A_solve_calls =
                        static_cast<std::size_t>(final_spectral.applications);
                }

                if (!run.converged)
                    ++algorithmic_nonconvergences;

                if (run.prediction_failure_residual)
                    ++prediction_failures_residual;
                if (run.prediction_failure_reference)
                    ++prediction_failures_reference;
                if (run.converged && !run.converged_to_reference)
                    ++wrong_equilibria;

                const std::size_t strategy_linear_calls =
                    run.picard_linear_solve_calls +
                    run.newton_linear_solve_calls +
                    run.spectral_A_solve_calls;

                const std::size_t diagnostic_initial_linear_calls =
                    (strategy == Strategy::PicardOnly ||
                     strategy == Strategy::NewtonOnly)
                        ? static_cast<std::size_t>(initial_spectral.applications)
                        : 0U;

                const std::size_t diagnostic_reference_linear_calls =
                    static_cast<std::size_t>(reference_spectral.applications);

                const std::size_t total_linear_calls_observed =
                    strategy_linear_calls +
                    diagnostic_initial_linear_calls +
                    diagnostic_reference_linear_calls +
                    diagnostic_spectral_A_solve_calls;

                const double seconds_case_total =
                    std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - case_wall_start).count();

                cases
                    << sample.case_id << ','
                    << sample.alpha << ','
                    << sample.theta_degrees << ','
                    << sample.rho_target << ','
                    << sample.rho_calibration_c1 << ','
                    << sample.c << ','
                    << method_name << ','
                    << mesh_name << ','
                    << actual_cells << ','
                    << strategy_name(strategy) << ','
                    << initial_guess_kind << ','
                    << run.rho_reference << ','
                    << run.rho_initial << ','
                    << run.rho_last_decision << ','
                    << run.rho_final_diagnostic << ','
                    << (run.rho_predicts_picard_convergence ? 1 : 0) << ','
                    << (run.converged ? 1 : 0) << ','
                    << (run.converged_to_reference ? 1 : 0) << ','
                    << (run.prediction_failure_residual ? 1 : 0) << ','
                    << (run.prediction_failure_reference ? 1 : 0) << ','
                    << (run.unexpected_picard_reference_convergence ? 1 : 0) << ','
                    << run.termination << ','
                    << run.outer_iterations << ','
                    << run.picard_steps << ','
                    << run.newton_steps << ','
                    << run.picard_linear_solve_calls << ','
                    << run.newton_linear_solve_calls << ','
                    << run.spectral_A_solve_calls << ','
                    << diagnostic_initial_linear_calls << ','
                    << diagnostic_reference_linear_calls << ','
                    << diagnostic_spectral_A_solve_calls << ','
                    << strategy_linear_calls << ','
                    << total_linear_calls_observed << ','
                    << run.krylov_iterations_total << ','
                    << run.newton_linear_hits_50 << ','
                    << run.rho_evaluations << ','
                    << run.solver_switches << ','
                    << run.eps_iterations_total << ','
                    << run.eps_applications_total << ','
                    << run.initial_residual << ','
                    << run.final_residual << ','
                    << run.initial_reference_error << ','
                    << run.final_reference_error << ','
                    << run.seconds_algorithmic << ','
                    << reference_rho_seconds << ','
                    << initial_rho_diagnostic_seconds << ','
                    << run.seconds_spectral_strategy << ','
                    << run.seconds_picard << ','
                    << run.seconds_newton << ','
                    << final_rho_seconds << ','
                    << case_setup_seconds << ','
                    << seconds_case_total << ','
                    << run.linear_iteration_limits << ','
                    << run.linear_residual_divergences << ','
                    << run.linear_other_failures << ','
                    << linear_reason_text(run.linear_reason_counts) << ','
                    << run.seconds_jacobian_assembly << ','
                    << run.regularisation_floor_active_iterations << '\n';
                cases.flush();

                metrics.push_back(CaseMetric{
                    std::string(method_name),
                    std::string(mesh_name),
                    run.rho_target,
                    run.rho_reference,
                    run.rho_initial,
                    run.converged,
                    false,  // case_exception
                    run.converged_to_reference,
                    run.prediction_failure_residual,
                    run.prediction_failure_reference,
                    run.unexpected_picard_reference_convergence,
                    strategy_linear_calls,
                    run.newton_linear_hits_50,
                    run.rho_evaluations,
                    run.solver_switches,
                    run.krylov_iterations_total,
                    run.seconds_algorithmic,
                    seconds_case_total,
                    initial_rho_diagnostic_seconds,
                    final_rho_seconds});

                if (settings.progress_every > 0U &&
                    sample.case_id % settings.progress_every == 0U) {
                    std::cout
                        << strategy_name(strategy) << " / "
                        << method_name << " / " << mesh_name
                        << " case " << sample.case_id
                        << "/" << settings.sample_count
                        << " rho_target=" << run.rho_target
                        << " rho_ref=" << run.rho_reference
                        << " rho0=" << run.rho_initial
                        << " conv=" << (run.converged ? 1 : 0)
                        << " ref=" << (run.converged_to_reference ? 1 : 0)
                        << " Ei=" << run.initial_reference_error
                        << " Ef=" << run.final_reference_error
                        << " NL=" << run.outer_iterations
                        << " LinCalls=" << strategy_linear_calls
                        << " t=" << run.seconds_algorithmic
                        << " s\n";
                }

                } catch (const std::exception& error) {
                    ++case_failures;

                    const double seconds_case_total =
                        std::chrono::duration<double>(
                            std::chrono::steady_clock::now()
                            - case_wall_start).count();

                    write_failed_case(
                        cases,
                        sample,
                        method_name,
                        mesh_name,
                        actual_cells,
                        strategy,
                        failure_initial_guess_kind,
                        error.what(),
                        seconds_case_total);

                    metrics.push_back(CaseMetric{
                        std::string(method_name),
                        std::string(mesh_name),
                        sample.rho_target,
                        std::numeric_limits<Real>::quiet_NaN(),
                        std::numeric_limits<Real>::quiet_NaN(),
                        false,   // converged
                        true,    // case_exception
                        false,   // converged_to_reference
                        false,   // prediction_failure_residual unavailable
                        false,   // prediction_failure_reference unavailable
                        false,   // unexpected Picard reference convergence
                        0U,      // linear calls
                        0U,      // linear hits 50
                        0U,      // rho evaluations
                        0U,      // solver switches
                        Index{0},
                        0.0,     // algorithmic seconds
                        seconds_case_total,
                        0.0,
                        0.0});

                    // This is deliberately not printed as SETUP/CAMPAIGN FAILURE.
                    // The loop proceeds to sample.case_id+1.
                    continue;
                }
            }

        } catch (const std::exception& error) {
            ++setup_failures;
            std::cerr
                << "\nSETUP/CAMPAIGN FAILURE: "
                << method_name << " / " << mesh_name
                << " : " << error.what() << '\n';
        }
    }
}

inline int run_program(
    Strategy strategy,
    const char* program_name) {

#ifndef GSMAKER_USE_SLEPC
    std::cerr
        << program_name
        << " requires GSMakerLib built with SLEPc\n";
    return 1;
#else
    try {
        gss::SlepcContext slepc;

        const auto settings =
            load_settings("config.yaml");
        auto random_cases =
            generate_random_cases(settings);

        const auto base_config =
            gss::YAMLConfig::load("config.yaml");
        auto study_settings =
            ps::load_study_settings("config.yaml");

        study_settings.random_mesh_seed = settings.mesh_seed;
        study_settings.requested_eigenpairs = settings.requested_eigenpairs;
        study_settings.eps_tolerance = settings.eps_tolerance;
        study_settings.eps_max_iterations = settings.eps_max_iterations;

        const auto result_dir =
            ps::create_result_directory();
        study::ConsoleTee console(
            result_dir / "console.log");

        std::ofstream samples(
            result_dir / "samples.csv");
        std::ofstream cases(
            result_dir / "cases.csv");
        std::ofstream switches(
            result_dir / "switches.csv");
        std::ofstream setup(
            result_dir / "setup.csv");
        std::ofstream manifest(
            result_dir / "manifest.txt");

        if (!samples || !cases || !switches || !setup || !manifest)
            throw std::runtime_error(
                "could not create result files");

        cases
            << std::scientific
            << std::setprecision(17)
            << "case_id,alpha,theta_degrees,rho_target,rho_calibration_c1,c,"
               "method,mesh_family,cells,strategy,initial_guess_kind,"
               "rho_reference,rho_initial,rho_last_decision,"
               "rho_final_diagnostic,rho_predicts_picard_convergence,"
               "converged,converged_to_reference,prediction_failure_residual,"
               "prediction_failure_reference,unexpected_picard_reference_convergence,"
               "termination,outer_iterations,picard_steps,newton_steps,"
               "picard_linear_solve_calls,newton_linear_solve_calls,"
               "strategy_spectral_A_solve_calls,diagnostic_initial_A_solve_calls,"
               "diagnostic_reference_A_solve_calls,diagnostic_final_A_solve_calls,"
               "strategy_linear_solve_calls,"
               "total_linear_solve_calls_observed,krylov_iterations_total,"
               "newton_linear_hits_50,rho_evaluations,solver_switches,"
               "eps_iterations_total,eps_applications_total,"
               "initial_residual,final_residual,initial_reference_error,"
               "final_reference_error,seconds_algorithmic,"
               "seconds_reference_rho_diagnostic,seconds_initial_rho_diagnostic,"
               "seconds_spectral_strategy,seconds_picard,seconds_newton,"
               "seconds_final_rho_diagnostic,seconds_case_setup,seconds_case_total,"
               "linear_iteration_limits,linear_residual_divergences,"
               "linear_other_failures,linear_reason_counts,seconds_jacobian_assembly,regularisation_floor_active_iterations\n";

        switches
            << std::scientific
            << std::setprecision(17)
            << "case_id,method,mesh_family,alpha,theta_degrees,rho_target,c,iteration,"
               "residual,rho,from_solver,to_solver\n";

        setup
            << std::scientific
            << std::setprecision(17)
            << "method,mesh_family,requested_cells,actual_cells,nominal_nx,"
               "nominal_ny,colours,jacobian_nnz,source_derivative,"
               "gradient_operator_relative_defect,calibration_rho_c1,"
               "spectral_A_factorisation_seconds,newton_A_factorisation_seconds\n";

        std::vector<CaseMetric> metrics;
        metrics.reserve(
            settings.sample_count *
            settings.mesh_families.size() * 2U);

        const vgt::gradients::polynomial::WLSQR2D wlsqr{
            vgt::gradients::WLSOptions2D{
                vgt::gradients::WLSWeighting2D::InverseDistanceSquared}};

        calibrate_case_scales(
            settings,
            random_cases,
            base_config,
            study_settings,
            wlsqr);

        const std::size_t calibration_failures =
            static_cast<std::size_t>(std::count_if(
                random_cases.begin(),
                random_cases.end(),
                [](const RandomCase& sample) {
                    return !sample.calibration_error.empty();
                }));

        samples << std::scientific << std::setprecision(17);
        samples
            << "case_id,alpha,theta_degrees,rho_target,rho_calibration_c1,c,"
               "calibration_error\n";
        for (const auto& sample : random_cases)
            samples
                << sample.case_id << ','
                << sample.alpha << ','
                << sample.theta_degrees << ','
                << sample.rho_target << ','
                << sample.rho_calibration_c1 << ','
                << sample.c << ','
                << csv_safe(sample.calibration_error) << '\n';
        samples.flush();

        const vgt::gradients::greengauss::ModifiedGreenGauss2D mgg{};

        std::size_t algorithmic_nonconvergences = 0U;
        std::size_t prediction_failures_residual = 0U;
        std::size_t prediction_failures_reference = 0U;
        std::size_t wrong_equilibria = 0U;
        std::size_t case_failures = 0U;
        std::size_t setup_failures = 0U;

        if(std::find(settings.methods.begin(),settings.methods.end(),"WLSQR")!=settings.methods.end())
        run_method(
            strategy,
            "WLSQR",
            wlsqr,
            settings,
            random_cases,
            base_config,
            study_settings,
            cases,
            switches,
            setup,
            metrics,
            algorithmic_nonconvergences,
            prediction_failures_residual,
            prediction_failures_reference,
            wrong_equilibria,
            case_failures,
            setup_failures);

        if(std::find(settings.methods.begin(),settings.methods.end(),"MGG")!=settings.methods.end())
        run_method(
            strategy,
            "MGG",
            mgg,
            settings,
            random_cases,
            base_config,
            study_settings,
            cases,
            switches,
            setup,
            metrics,
            algorithmic_nonconvergences,
            prediction_failures_residual,
            prediction_failures_reference,
            wrong_equilibria,
            case_failures,
            setup_failures);

        write_summary(result_dir, metrics, settings);
        write_rho_bins(result_dir, metrics, settings.bins_by_target);

        manifest
            << "program=" << program_name << '\n'
            << "paper=PaperIII\n"
            << "campaign=random_four_solver_strategies\n"
            << "strategy=" << strategy_name(strategy) << '\n'
            << "sample_count=" << settings.sample_count << '\n'
            << "sampling_seed=" << settings.sampling_seed << '\n'
            << "sampling_generator=std::mt19937_64\n"
            << "sampling_population_generated_before_any_nonlinear_solve=1\n"
            << "sampling_stratification=case_id_modulo_5\n"
            << "rho_target_strata=[0.70,0.90);[0.90,0.98);[0.98,1.02);"
               "[1.02,1.10);[1.10,1.30]\n"
            << "rho_target_cases_per_stratum="
            << (settings.sample_count + 4U) / 5U << ','
            << (settings.sample_count + 3U) / 5U << ','
            << (settings.sample_count + 2U) / 5U << ','
            << (settings.sample_count + 1U) / 5U << ','
            << settings.sample_count / 5U << '\n' 
            << "rho_target_calibration=Cartesian_WLSQR_reference_state_c1_scaling\n"
            << "rho_target_calibration_failures=" << calibration_failures << '\n'
            << "alpha_range=(" << settings.alpha_min << ',' << settings.alpha_max << ")\n"
            << "theta_range_degrees=[" << settings.theta_min_degrees << ','
            << settings.theta_max_degrees << "]\n"
            << "requested_cells=" << settings.target_cells << '\n'
            << "nominal_nx=" << settings.nominal_nx << '\n'
            << "nominal_ny=" << settings.nominal_ny << '\n'
            << "mesh_families=" << [&] {std::string t; for(const auto& f:settings.mesh_families) {if(!t.empty())t+=","; t+=f;} return t;}() << '\n' 
            << "methods=" << [&] {std::string t; for(const auto& f:settings.methods) {if(!t.empty())t+=","; t+=f;} return t;}() << '\n'
            << "diagnostics.reference_rho=" << settings.reference_rho << '\n'
            << "diagnostics.initial_rho=" << settings.initial_rho << '\n'
            << "diagnostics.jacobian_assembly=" << settings.jacobian_timing << '\n'
            << "diagnostics.regularisation_floor=" << settings.floor_activity << '\n'
            << "rho_bins_basis=" << (settings.bins_by_target?"rho_target":"rho_reference") << '\n'
            << "seconds_jacobian_assembly_scope=D_assembly_consumed_by_Newton_excluding_Picard_spectral_diagnostics\n"
            << "regularisation_floor_definition=at_least_one_cell_g2_le_delta_squared_times_max_g2_at_attempted_nonlinear_step\n"
            << "initial_rho_if_disabled=mixed_strategy_still_computes_required_selector_rho_but_column_is_nan\n"
            << "serial_optimisation=invariant_workspace_cached_source_inplace_J\n"
            << "newton_initial_B_reuse=assembly_time_charged_to_algorithm\n"
            << "sparse_storage=adjacency_CSR\n"
            << "adjacency_storage_runtime_gate=enabled\n"
            << "source_derivative=" << settings.source_derivative << '\n'
            << "source_jacobian="
            << (settings.source_derivative == "adolc"
                    ? "adolc_coloured_forward_mode"
                    : "coloured_sparse_finite_difference_explicit")
            << '\n'
            << "finite_difference_fallback=disabled\n"
            << "gradient_operator_extraction=exact_coloured_superposition\n"
            << "regularisation_threshold_linearisation=frozen_at_current_state\n"
            << "spectral_operator=matrix_free_A_inverse_B\n"
            << "newton_linear_solver=GMRES\n"
            << "newton_preconditioner=invariant_A_factorisation\n"
            << "linear_relative_tolerance=" << settings.krylov_relative_tolerance << '\n'
            << "linear_max_iterations=" << settings.krylov_max_iterations << '\n'
            << "linear_nonconvergence=nonfatal_record_reason_use_best_iterate\n"
            << "linear_iteration_limit_reason=KSP_DIVERGED_ITS_not_counted_as_divergence\n"
            << "linear_residual_divergence_reason=KSP_DIVERGED_DTOL\n"
            << "linear_reason_counts_format=code:name:count_separated_by_semicolons\n"
            << "linear_status_scope=completed_Newton_KSPSolve_calls_per_case\n"
            << "nonlinear_relative_tolerance=" << settings.tolerance_relative << '\n'
            << "nonlinear_absolute_tolerance=" << settings.tolerance_absolute << '\n'
            << "nonlinear_max_iterations=" << settings.max_iterations << '\n'
            << "initialisation=reference_plus_deterministic_perturbation\n"
            << "perturbation_amplitude=" << settings.perturbation_amplitude << '\n'
            << "continuation=disabled_independent_random_cases\n"
            << "reference_relative_tolerance="
            << settings.reference_relative_tolerance << '\n'
            << "algorithmic_nonconvergences="
            << algorithmic_nonconvergences << '\n'
            << "prediction_failures_residual="
            << prediction_failures_residual << '\n'
            << "prediction_failures_reference="
            << prediction_failures_reference << '\n'
            << "converged_wrong_equilibria="
            << wrong_equilibria << '\n'
            << "case_exceptions=" << case_failures << '\n'
            << "setup_failures=" << setup_failures << '\n'
            << "outputs=samples.csv;cases.csv;switches.csv;setup.csv;"
               "summary.csv;rho_bins.csv;console.log;manifest.txt\n";

        std::cout
            << "\nResults: " << result_dir.string() << '\n'
            << "Cases recorded: " << metrics.size() << '\n'
            << "rho-target calibration failures: " << calibration_failures << '\n'
            << "Algorithmic non-convergences: "
            << algorithmic_nonconvergences << '\n'
            << "rho_reference<1 but residual did not converge: "
            << prediction_failures_residual << '\n'
            << "rho_reference<1 but reference equilibrium not recovered: "
            << prediction_failures_reference << '\n'
            << "Converged to a different equilibrium: "
            << wrong_equilibria << '\n'
            << "Case exceptions (campaign continued): "
            << case_failures << '\n'
            << "Setup failures: " << setup_failures << '\n';

        return setup_failures == 0U ? 0 : 1;

    } catch (const std::exception& error) {
        std::cerr
            << program_name
            << " failed: "
            << error.what() << '\n';
        return 1;
    }
#endif
}

} // namespace p3official
