#pragma once

#include <adolc/adolc.h>
#include <adolc/adolc_sparse.h>

#include <cstddef>
#include <cstdlib>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gss::autodiff {

struct SparseJacobianEntry {
    std::size_t row{};
    std::size_t col{};
    double value{};
};

struct SparseJacobianResult {
    std::size_t rows{};
    std::size_t cols{};
    std::vector<SparseJacobianEntry> entries;
};

class AdolcSparseJacobian {
public:
    explicit AdolcSparseJacobian(short tag) : tag_(tag) {}

    [[nodiscard]] short tag() const noexcept { return tag_; }

    template <class ActiveFunction>
    void trace(std::span<const double> x0,
               std::size_t output_size,
               ActiveFunction&& function) const {

        const int n = checkedInt(x0.size());
        const int m = checkedInt(output_size);

        trace_on(tag_);

        std::vector<adouble> x(static_cast<std::size_t>(n));
        std::vector<adouble> y(static_cast<std::size_t>(m));
        std::vector<double> y0(static_cast<std::size_t>(m), 0.0);

        for (int i=0;i<n;++i) {
            x[static_cast<std::size_t>(i)] <<= x0[static_cast<std::size_t>(i)];
        }

        function(std::span<const adouble>(x.data(), x.size()),
                 std::span<adouble>(y.data(), y.size()));

        for (int i=0;i<m;++i) {
            y[static_cast<std::size_t>(i)] >>= y0[static_cast<std::size_t>(i)];
        }

        trace_off();
    }

    [[nodiscard]] SparseJacobianResult evaluate(
        std::span<const double> x,
        std::size_t output_size,
        bool repeat = false) const {

        const int n = checkedInt(x.size());
        const int m = checkedInt(output_size);

        std::vector<double> point(x.begin(), x.end());
        int nnz = 0;
        unsigned int* rind = nullptr;
        unsigned int* cind = nullptr;
        double* values = nullptr;

        // ADOL-C sparse_jac options:
        // options[0] = 0, index-domain sparsity identification
        // options[1] = 0, safe mode
        // options[2] = 0, automatic compression choice
        // options[3] = 0, column compression when selected
        int options[4] = {0,0,0,0};

        const int rc = sparse_jac(
            tag_, m, n, repeat ? 1 : 0, point.data(),
            &nnz, &rind, &cind, &values, options);

        if (rc < 0) {
            std::free(rind);
            std::free(cind);
            std::free(values);
            throw std::runtime_error("ADOL-C sparse_jac failed");
        }

        SparseJacobianResult result;
        result.rows = static_cast<std::size_t>(m);
        result.cols = static_cast<std::size_t>(n);
        result.entries.reserve(static_cast<std::size_t>(nnz));

        for (int k=0;k<nnz;++k) {
            result.entries.push_back({
                static_cast<std::size_t>(rind[k]),
                static_cast<std::size_t>(cind[k]),
                values[k]});
        }

        std::free(rind);
        std::free(cind);
        std::free(values);
        return result;
    }

private:
    static int checkedInt(std::size_t n) {
        if (n > static_cast<std::size_t>(2147483647)) {
            throw std::overflow_error("ADOL-C dimension exceeds int range");
        }
        return static_cast<int>(n);
    }

    short tag_{};
};

} // namespace gss::autodiff
