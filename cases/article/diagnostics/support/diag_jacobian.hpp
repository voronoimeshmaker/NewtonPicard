#pragma once
// Copied from support/official.hpp lines 457-830, commit b8b5f057d68d5ecdcbf9c4a3abe9f6155a26854a.
// Only class name and access changed; production arithmetic is unchanged.
namespace p3official {
template<class GradientMethod>
class DiagJacobian final {
public:
    DiagJacobian(
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
                "DiagJacobian: unsupported derivative method");

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

    DiagJacobian(const DiagJacobian&) = delete;
    DiagJacobian& operator=(const DiagJacobian&) = delete;

    template<class Source>
    [[nodiscard]] std::vector<Real> values(
        std::span<const Real> state,
        const Source& source) {

        if (state.size() != pattern_.cells)
            throw std::invalid_argument(
                "DiagJacobian: state size does not match mesh");

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

public:
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

}
