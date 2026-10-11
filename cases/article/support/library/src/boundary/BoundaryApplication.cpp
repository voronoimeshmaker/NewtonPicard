#include "gsmaker/boundary/BoundaryApplication.hpp"

namespace gss {

std::vector<EvaluatedBoundaryFace> BoundaryApplication::evaluate(
    std::span<const BoundaryFaceGeometry> faces,
    const BoundaryCondition& condition) {

    std::vector<EvaluatedBoundaryFace> result;
    result.reserve(faces.size());

    for (const auto& face : faces) {
        const BoundaryPoint2D point{
            face.owner_cell,
            face.x,
            face.y,
            face.normal_x,
            face.normal_y};

        result.push_back({face, coefficients_at(condition, point)});
    }

    return result;
}

} // namespace gss

