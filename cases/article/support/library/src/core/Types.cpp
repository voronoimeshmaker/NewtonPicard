#include "gsmaker/core/Types.hpp"

namespace gss {

void CellData::resize(Index n) {
    n_cells = n;

    const auto s =
        static_cast<std::size_t>(n);

    psi.assign(s, Real{0});
    psi_old.assign(s, Real{0});
    source.assign(s, Real{0});

    grad_psi_x.assign(s, Real{0});
    grad_psi_y.assign(s, Real{0});

    grad_F_x.assign(s, Real{0});
    grad_F_y.assign(s, Real{0});

    kappa.assign(s, Real{1});
}

} // namespace gss
