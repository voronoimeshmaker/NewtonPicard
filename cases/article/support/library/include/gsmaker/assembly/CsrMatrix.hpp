#pragma once

#include <span>
#include <vector>

#include "gsmaker/core/Types.hpp"



namespace gss {

/// @brief Positive mathematical properties known for a CSR matrix.
///
/// A false value means "not asserted", not "known to be false".  This keeps
/// CsrMatrix generic while allowing assemblers that know more about their
/// operator to communicate that information to downstream linear algebra.
struct CsrMatrixProperties {
    bool structurally_symmetric{false};
    bool symmetric{false};
    bool positive_definite{false};
};

class CsrMatrix final {
public:
    CsrMatrix() = default;
    CsrMatrix(
        Index rows,
        std::vector<Index> row_start,
        std::vector<Index> column,
        std::vector<Real> value,
        CsrMatrixProperties properties = {});

    [[nodiscard]] Index rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t nnz() const noexcept { return value_.size(); }

    [[nodiscard]] std::span<const Index> row_start() const noexcept { return row_start_; }
    [[nodiscard]] std::span<const Index> column() const noexcept { return column_; }
    [[nodiscard]] std::span<const Real> value() const noexcept { return value_; }
    [[nodiscard]] const CsrMatrixProperties& properties() const noexcept { return properties_; }

    [[nodiscard]] std::vector<Real> multiply(std::span<const Real> x) const;
    [[nodiscard]] Real coefficient(Index row, Index column) const;
    void validate() const;

private:
    Index rows_{0};
    std::vector<Index> row_start_{};
    std::vector<Index> column_{};
    std::vector<Real> value_{};
    CsrMatrixProperties properties_{};
};

} // namespace gss


