#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include "gsmaker/core/MeshAdapter.hpp"
#include "gsmaker/core/Types.hpp"
#include "gsmaker/semilinear/GradientRatioSource.hpp"

namespace gss {

enum class RecoveryMode {
    None,
    WeightedAverage,
    LinearInPsi,
    QuadraticInPsi
};

struct StationaryPointRecoveryResult {
    /// Source derivative with the flagged cells repaired.
    std::vector<Real> values;

    /// Mode actually used at each cell; None where nothing was done.
    std::vector<RecoveryMode> mode;

    std::size_t recovered_cells{};
    std::size_t unrecoverable_cells{};
    std::size_t sweeps{};

    /// How many repairs used each mode.
    std::size_t averaged{};
    std::size_t linear_fits{};
    std::size_t quadratic_fits{};
};

/// @brief Repair the source derivative at stationary points of Psi.
///
/// Where grad Psi vanishes the gradient ratio is 0/0 and carries no
/// information, however the regularisation threshold is chosen.  The value
/// there is not small: for the Fadeev solution the stationary point at the
/// origin has Psi = log(1.25) and dF/dPsi = 0.64, so a flagged cell that keeps
/// its regularised value is wrong by a hundred per cent at a single point.
/// That error enters the row and column of the source Jacobian belonging to
/// that cell, and therefore the spectrum of the Picard iteration operator.
///
/// THE REPAIR EXPLOITS THAT dF/dPsi IS A FUNCTION OF Psi.  Every reliable
/// neighbour supplies a pair (Psi_k, r_k) on the same curve r = F'(Psi).  The
/// repaired value is obtained by fitting that curve over the neighbours, in
/// the shifted variable s = Psi_k - Psi_P, and taking the constant term, which
/// is the fit evaluated at Psi_P.  Fitting in the shifted variable means the
/// answer is read off directly rather than extrapolated.
///
/// WHY NOT SIMPLY AVERAGE.  An average is the degree-zero case and is only
/// correct when Psi is nearly constant over the stencil.  That holds on a fine
/// mesh, where Psi_k - Psi_P = O(h^2) at a stationary point, and fails badly on
/// a coarse one: on the 211-cell hexagonal mesh the neighbours of the origin
/// carry dF/dPsi between 0.23 and 0.81 against a true value of 0.64, and their
/// average is 0.38.  The fit covers both regimes, and falls back to the average
/// exactly when the spread in Psi is too small for a fit to be conditioned,
/// which is the regime in which the average is the right answer.
///
/// The Fadeev stationary point is a saddle, not an extremum, so Psi_P lies
/// inside the range spanned by the neighbours and the fit interpolates.  The
/// implementation does not rely on that: the fitted value is confined to the
/// range of the neighbour values, so an ill-posed extrapolation cannot produce
/// a value the data does not support.
///
/// Cells adjacent to one another may both be flagged on a fine mesh, so the
/// repair is swept: a cell recovered in one sweep becomes a reliable donor in
/// the next.
class StationaryPointRecovery final {
public:
    [[nodiscard]]
    static StationaryPointRecoveryResult apply(
        const MeshTopology& topology,
        const GradientRatioSourceResult& source,
        std::span<const Real> psi,
        std::size_t max_sweeps = 5U) {

        const std::size_t n = source.values.size();

        if (n == 0U || source.regularised.size() != n) {
            throw std::invalid_argument(
                "StationaryPointRecovery requires a source result with one flag per cell");
        }
        if (psi.size() != n) {
            throw std::invalid_argument(
                "StationaryPointRecovery psi must have one value per cell");
        }
        if (topology.coord_x.size() != n || topology.coord_y.size() != n) {
            throw std::invalid_argument(
                "StationaryPointRecovery topology and source sizes must agree");
        }
        if (max_sweeps == 0U) {
            throw std::invalid_argument(
                "StationaryPointRecovery requires at least one sweep");
        }

        StationaryPointRecoveryResult result;
        result.values = source.values;
        result.mode.assign(n, RecoveryMode::None);

        std::vector<std::uint8_t> pending = source.regularised;
        std::vector<std::uint8_t> reliable(n, std::uint8_t{0});

        for (std::size_t i = 0; i < n; ++i) {
            const bool flagged = pending[i] != std::uint8_t{0};
            reliable[i] = flagged ? std::uint8_t{0} : std::uint8_t{1};
        }

        for (std::size_t sweep = 0; sweep < max_sweeps; ++sweep) {

            std::vector<std::size_t> index;
            std::vector<Real> value;
            std::vector<RecoveryMode> used;

            for (std::size_t i = 0; i < n; ++i) {

                if (pending[i] == std::uint8_t{0}) {
                    continue;
                }

                const Fit fit = fit_cell(topology, result.values, psi, reliable, i, n);

                if (fit.mode == RecoveryMode::None) {
                    continue;
                }

                index.push_back(i);
                value.push_back(fit.value);
                used.push_back(fit.mode);
            }

            ++result.sweeps;

            if (index.empty()) {
                break;
            }

            // Applied only after every flagged cell has been visited, so the
            // outcome cannot depend on the order of traversal.
            for (std::size_t j = 0; j < index.size(); ++j) {

                const std::size_t i = index[j];

                result.values[i] = value[j];
                result.mode[i] = used[j];
                reliable[i] = std::uint8_t{1};
                pending[i] = std::uint8_t{0};

                ++result.recovered_cells;

                if (used[j] == RecoveryMode::QuadraticInPsi) {
                    ++result.quadratic_fits;
                } else if (used[j] == RecoveryMode::LinearInPsi) {
                    ++result.linear_fits;
                } else {
                    ++result.averaged;
                }
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            if (pending[i] != std::uint8_t{0}) {
                ++result.unrecoverable_cells;
            }
        }

        return result;
    }

private:
    struct Fit {
        Real value{};
        RecoveryMode mode{RecoveryMode::None};
    };

    struct Sample {
        double s{};
        double r{};
        double w{};
    };

    /// Relative spread in Psi below which a fit is not attempted and the
    /// weighted average is used instead.  At that point Psi is flat over the
    /// stencil and the average is the correct degree-zero answer.
    static constexpr double minimum_relative_spread = 1.0e-6;

    [[nodiscard]]
    static Fit fit_cell(
        const MeshTopology& topology,
        const std::vector<Real>& values,
        std::span<const Real> psi,
        const std::vector<std::uint8_t>& reliable,
        std::size_t cell,
        std::size_t n) {

        const auto begin = static_cast<std::size_t>(topology.neighbor_start[cell]);
        const auto end = static_cast<std::size_t>(topology.neighbor_start[cell + 1U]);

        const double xi = static_cast<double>(topology.coord_x[cell]);
        const double yi = static_cast<double>(topology.coord_y[cell]);
        const double psi_i = static_cast<double>(psi[cell]);

        std::vector<Sample> sample;
        sample.reserve(end - begin);

        double lowest = 0.0;
        double highest = 0.0;
        double psi_lowest = 0.0;
        double psi_highest = 0.0;
        bool first = true;

        for (std::size_t k = begin; k < end; ++k) {

            const auto neighbour = static_cast<std::size_t>(topology.neighbor_list[k]);

            if (neighbour >= n) {
                continue;
            }
            if (reliable[neighbour] == std::uint8_t{0}) {
                continue;
            }

            const double dx = static_cast<double>(topology.coord_x[neighbour]) - xi;
            const double dy = static_cast<double>(topology.coord_y[neighbour]) - yi;
            const double distance_squared = dx * dx + dy * dy;

            if (!(distance_squared > 0.0)) {
                continue;
            }

            Sample entry;
            entry.s = static_cast<double>(psi[neighbour]) - psi_i;
            entry.r = static_cast<double>(values[neighbour]);
            entry.w = 1.0 / distance_squared;

            if (!std::isfinite(entry.s) || !std::isfinite(entry.r)) {
                continue;
            }

            if (first) {
                lowest = entry.r;
                highest = entry.r;
                psi_lowest = entry.s;
                psi_highest = entry.s;
                first = false;
            } else {
                lowest = std::min(lowest, entry.r);
                highest = std::max(highest, entry.r);
                psi_lowest = std::min(psi_lowest, entry.s);
                psi_highest = std::max(psi_highest, entry.s);
            }

            sample.push_back(entry);
        }

        Fit fit;

        if (sample.empty()) {
            return fit;
        }

        const double spread = psi_highest - psi_lowest;
        const double psi_scale = std::max(std::abs(psi_i), 1.0);
        const bool spread_is_usable = (spread / psi_scale) > minimum_relative_spread;

        if (sample.size() >= 3U && spread_is_usable) {
            const Fit quadratic = solve_polynomial(sample, spread, 2U);
            if (quadratic.mode != RecoveryMode::None) {
                return clamp_to_data(quadratic, lowest, highest, RecoveryMode::QuadraticInPsi);
            }
        }

        if (sample.size() >= 2U && spread_is_usable) {
            const Fit linear = solve_polynomial(sample, spread, 1U);
            if (linear.mode != RecoveryMode::None) {
                return clamp_to_data(linear, lowest, highest, RecoveryMode::LinearInPsi);
            }
        }

        double weighted = 0.0;
        double total = 0.0;

        for (const Sample& entry : sample) {
            weighted += entry.w * entry.r;
            total += entry.w;
        }

        if (!(total > 0.0)) {
            return fit;
        }

        const double average = weighted / total;

        if (!std::isfinite(average)) {
            return fit;
        }

        fit.value = static_cast<Real>(average);
        fit.mode = RecoveryMode::WeightedAverage;
        return fit;
    }

    /// Weighted least squares of r against powers of the scaled shift
    /// t = s / spread, returning the constant term, which is the fit at
    /// t = 0, that is at Psi_P.  Scaling by the spread keeps the normal
    /// matrix well conditioned whatever the magnitude of Psi.
    [[nodiscard]]
    static Fit solve_polynomial(
        const std::vector<Sample>& sample,
        double spread,
        std::size_t degree) {

        Fit fit;

        if (!(spread > 0.0)) {
            return fit;
        }

        const std::size_t size = degree + 1U;

        double matrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
        double rhs[3] = {0.0, 0.0, 0.0};

        for (const Sample& entry : sample) {

            const double t = entry.s / spread;
            double basis[3] = {1.0, t, t * t};

            for (std::size_t a = 0; a < size; ++a) {
                for (std::size_t b = 0; b < size; ++b) {
                    matrix[a][b] += entry.w * basis[a] * basis[b];
                }
                rhs[a] += entry.w * basis[a] * entry.r;
            }
        }

        double solution[3] = {0.0, 0.0, 0.0};

        if (!solve_symmetric(matrix, rhs, solution, size)) {
            return fit;
        }

        if (!std::isfinite(solution[0])) {
            return fit;
        }

        fit.value = static_cast<Real>(solution[0]);
        fit.mode = RecoveryMode::WeightedAverage;
        return fit;
    }

    /// Gaussian elimination with partial pivoting on a small dense system.
    /// Returns false when the system is too close to singular to trust.
    [[nodiscard]]
    static bool solve_symmetric(
        double matrix[3][3],
        const double rhs[3],
        double solution[3],
        std::size_t size) {

        double a[3][4] = {{0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};

        double scale = 0.0;

        for (std::size_t i = 0; i < size; ++i) {
            for (std::size_t j = 0; j < size; ++j) {
                a[i][j] = matrix[i][j];
                scale = std::max(scale, std::abs(matrix[i][j]));
            }
            a[i][size] = rhs[i];
        }

        if (!(scale > 0.0)) {
            return false;
        }

        const double tolerance = 1.0e-12 * scale;

        for (std::size_t column = 0; column < size; ++column) {

            std::size_t pivot = column;

            for (std::size_t row = column + 1U; row < size; ++row) {
                if (std::abs(a[row][column]) > std::abs(a[pivot][column])) {
                    pivot = row;
                }
            }

            if (std::abs(a[pivot][column]) < tolerance) {
                return false;
            }

            if (pivot != column) {
                for (std::size_t j = 0; j <= size; ++j) {
                    const double swap = a[column][j];
                    a[column][j] = a[pivot][j];
                    a[pivot][j] = swap;
                }
            }

            for (std::size_t row = column + 1U; row < size; ++row) {
                const double factor = a[row][column] / a[column][column];
                for (std::size_t j = column; j <= size; ++j) {
                    a[row][j] -= factor * a[column][j];
                }
            }
        }

        for (std::size_t k = size; k-- > 0;) {
            double sum = a[k][size];
            for (std::size_t j = k + 1U; j < size; ++j) {
                sum -= a[k][j] * solution[j];
            }
            solution[k] = sum / a[k][k];
            if (!std::isfinite(solution[k])) {
                return false;
            }
        }

        return true;
    }

    /// A fit may not produce a value the neighbours do not support.  F' is
    /// smooth, so the repaired value belongs inside the range the stencil
    /// spans; a fit that leaves it has been driven by conditioning rather than
    /// by data.
    [[nodiscard]]
    static Fit clamp_to_data(Fit fit, double lowest, double highest, RecoveryMode mode) {
        const double raw = static_cast<double>(fit.value);
        const double bounded = std::min(std::max(raw, lowest), highest);
        fit.value = static_cast<Real>(bounded);
        fit.mode = mode;
        return fit;
    }
};

} // namespace gss