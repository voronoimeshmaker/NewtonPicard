#pragma once

#include <span>
#include <vector>


#include "gsmaker/assembly/CsrMatrix.hpp"
#include "gsmaker/boundary/BoundaryApplication.hpp"
#include "gsmaker/core/MeshAdapter.hpp"

namespace gss {

struct PoissonSystem {
    CsrMatrix matrix;
    std::vector<Real> boundary_rhs;
};

/// @brief Assemble the integrated finite-volume operator for -Laplacian(Psi).
class PoissonAssembler final {
public:
    [[nodiscard]] static PoissonSystem assemble(
        const MeshAdapter& mesh,
        std::span<const EvaluatedBoundaryFace> boundary,
        Real diffusivity = Real{1});

    /// @brief Build A Psi = -V source + boundary_rhs for Laplacian(Psi)=source.
    [[nodiscard]] static std::vector<Real> build_rhs(
        const MeshAdapter& mesh,
        std::span<const Real> source,
        std::span<const EvaluatedBoundaryFace> boundary,
        Real diffusivity = Real{1});
};

} // namespace gss


