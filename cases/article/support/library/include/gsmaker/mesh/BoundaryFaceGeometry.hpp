#pragma once

#include "gsmaker/core/Types.hpp"

namespace gss {

/// @brief Solver-independent geometry of one physical boundary face.
struct BoundaryFaceGeometry {
    Index owner_cell{};
    Real x{};
    Real y{};
    Real normal_x{};
    Real normal_y{};
    Real measure{};
    Real perpendicular_distance{};
    bool application_point_inside_face{true};
};

} // namespace gss


