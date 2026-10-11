#include "gsmaker/assembly/CsrMatrix.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>




namespace gss {

CsrMatrix::CsrMatrix(
    Index rows,
    std::vector<Index> row_start,
    std::vector<Index> column,
    std::vector<Real> value,
    CsrMatrixProperties properties)
    : rows_(rows),
      row_start_(std::move(row_start)),
      column_(std::move(column)),
      value_(std::move(value)),
      properties_(properties) {

    // SPD is, by definition here, a stronger positive assertion than
    // symmetry; numerical symmetry in turn implies structural symmetry.
    if (properties_.positive_definite) {
        properties_.symmetric = true;
    }
    if (properties_.symmetric) {
        properties_.structurally_symmetric = true;
    }

    validate();
}

void CsrMatrix::validate() const {
    if (rows_ < 0) throw std::runtime_error("CsrMatrix row count must be non-negative");
    const auto rows = static_cast<std::size_t>(rows_);
    if (row_start_.size() != rows + 1U) {
        throw std::runtime_error("CsrMatrix row_start must contain rows + 1 entries");
    }
    if (column_.size() != value_.size()) {
        throw std::runtime_error("CsrMatrix column and value arrays must have equal size");
    }
    if (row_start_.empty() || row_start_.front() != Index{0} ||
        row_start_.back() != static_cast<Index>(value_.size())) {
        throw std::runtime_error("CsrMatrix row offsets are inconsistent");
    }
    for (std::size_t r = 0; r < rows; ++r) {
        if (row_start_[r] > row_start_[r + 1U]) {
            throw std::runtime_error("CsrMatrix row offsets must be monotone");
        }
    }
    for (std::size_t k = 0; k < column_.size(); ++k) {
        if (column_[k] < 0 || column_[k] >= rows_) {
            throw std::runtime_error("CsrMatrix contains an out-of-range column");
        }
        if (!std::isfinite(static_cast<double>(value_[k]))) {
            throw std::runtime_error("CsrMatrix contains a non-finite coefficient");
        }
    }
}

std::vector<Real> CsrMatrix::multiply(std::span<const Real> x) const {
    if (x.size() != static_cast<std::size_t>(rows_)) {
        throw std::invalid_argument("CsrMatrix multiply vector size mismatch");
    }
    std::vector<Real> result(x.size(), Real{0});
    for (Index r = 0; r < rows_; ++r) {
        Real sum = Real{0};
        const auto rs = static_cast<std::size_t>(r);
        for (Index k = row_start_[rs]; k < row_start_[rs + 1U]; ++k) {
            const auto ks = static_cast<std::size_t>(k);
            sum += value_[ks] * x[static_cast<std::size_t>(column_[ks])];
        }
        result[rs] = sum;
    }
    return result;
}

Real CsrMatrix::coefficient(Index row, Index column) const {
    if (row < 0 || row >= rows_ || column < 0 || column >= rows_) {
        throw std::out_of_range("CsrMatrix coefficient index out of range");
    }
    const auto rs = static_cast<std::size_t>(row);
    for (Index k = row_start_[rs]; k < row_start_[rs + 1U]; ++k) {
        const auto ks = static_cast<std::size_t>(k);
        if (column_[ks] == column) return value_[ks];
    }
    return Real{0};
}

} // namespace gss
