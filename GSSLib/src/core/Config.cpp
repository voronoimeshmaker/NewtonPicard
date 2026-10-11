#include "gss/core/Types.hpp"

namespace gss {

void CellData::resize(Index n) {
    n_cells = n;
    const auto s = static_cast<std::size_t>(n);
    psi.assign(s, 0.0);
    psi_old.assign(s, 0.0);
    source.assign(s, 0.0);
    volume.assign(s, 0.0);
    grad_psi_x.assign(s, 0.0);
    grad_psi_y.assign(s, 0.0);
    grad_F_x.assign(s, 0.0);
    grad_F_y.assign(s, 0.0);
    kappa.assign(s, 1.0);
}

} // namespace gss
