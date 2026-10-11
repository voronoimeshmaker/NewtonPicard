#include "gsmaker/solver/SolutionState.hpp"

#include <stdexcept>

namespace gss {

namespace {

[[nodiscard]] std::size_t checked_size(Index n_cells) {
    if (n_cells < Index{0}) {
        throw std::invalid_argument(
            "SolutionState cell count must be non-negative");
    }

    return static_cast<std::size_t>(n_cells);
}

} // namespace


SolutionState::SolutionState(Index n) {
    resize(n);
}



void SolutionState::resize(Index n) {
    const auto size = checked_size(n);

    n_cells = n;

    psi.assign(size, Real{0});
    psi_old.assign(size, Real{0});
    source.assign(size, Real{0});

    grad_psi_x.assign(size, Real{0});
    grad_psi_y.assign(size, Real{0});
    grad_F_x.assign(size, Real{0});
    grad_F_y.assign(size, Real{0});

    kappa.assign(size, Real{1});
}


void SolutionState::validate() const {
    const auto expected = checked_size(n_cells);

    const auto require_size =
        [expected](
            const std::vector<Real>& field,
            const char* name) {

            if (field.size() != expected) {
                throw std::runtime_error(
                    std::string("SolutionState field '") +
                    name +
                    "' does not match cell count");
            }
        };

    require_size(psi, "psi");
    require_size(psi_old, "psi_old");
    require_size(source, "source");

    require_size(grad_psi_x, "grad_psi_x");
    require_size(grad_psi_y, "grad_psi_y");
    require_size(grad_F_x, "grad_F_x");
    require_size(grad_F_y, "grad_F_y");

    require_size(kappa, "kappa");
}

} // namespace gss
