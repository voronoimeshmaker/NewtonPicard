#pragma once

#include <span>
#include <vector>

#include "gsmaker/boundary/BoundaryCondition.hpp"
#include "gsmaker/mesh/BoundaryFaceGeometry.hpp"


namespace gss {

struct EvaluatedBoundaryFace {
    BoundaryFaceGeometry geometry{};
    BoundaryCoefficients coefficients{};
};

class BoundaryApplication final {
public:
    [[nodiscard]] static std::vector<EvaluatedBoundaryFace> evaluate(
        std::span<const BoundaryFaceGeometry> faces,
        const BoundaryCondition& condition);
};

} // namespace gss
