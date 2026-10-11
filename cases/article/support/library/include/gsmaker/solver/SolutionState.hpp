#pragma once

#include <vector>

#include "gsmaker/core/Types.hpp"

namespace gss {

/**
 * @brief Cell-centred numerical state used by the nonlinear solver.
 *
 * This type contains only solution-dependent or solution-derived quantities.
 * It intentionally contains no mesh geometry.
 *
 * The member names and initial values mirror the corresponding legacy
 * CellData fields to minimise migration risk.
 */
struct SolutionState {
    Index n_cells{0};

    std::vector<Real> psi{};
    std::vector<Real> psi_old{};
    std::vector<Real> source{};

    std::vector<Real> grad_psi_x{};
    std::vector<Real> grad_psi_y{};
    std::vector<Real> grad_F_x{};
    std::vector<Real> grad_F_y{};

    std::vector<Real> kappa{};


    
    SolutionState() = default;
    explicit SolutionState(Index n);

    /**
     * @brief Reset all state fields for n cells.
     *
     * This preserves the behaviour of the legacy CellData::resize:
     * solution/source/gradient fields are reset to zero and kappa to one.
     */
    void resize(Index n);

    /**
     * @brief Check that every state array has exactly n_cells entries.
     */
    void validate() const;
};

} // namespace gss
