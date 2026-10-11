#pragma once

#include <stdexcept>
#include <vector>

#include <petsc.h>

#include "gss/core/Config.hpp"
#include "gss/core/Types.hpp"

namespace gss {

/// @brief RAII helper for PETSc vectors and matrix-vector operations.
class PETScLinearSystem {
public:
    /// @brief Create a distributed PETSc vector initialized from std::vector data.
    [[nodiscard]] static Vec make_vec(std::span<const Real> values, MPI_Comm comm = PETSC_COMM_WORLD);

    /// @brief Create a distributed PETSc vector with the requested global length.
    [[nodiscard]] static Vec make_vec(Index global_size, MPI_Comm comm = PETSC_COMM_WORLD);

    /// @brief Copy PETSc vector values to std::vector.
    [[nodiscard]] static std::vector<Real> to_std_vector(Vec vec);

    /// @brief Destroy a PETSc vector if non-null.
    static void destroy(Vec& vec) noexcept;

    /// @brief Destroy a PETSc matrix if non-null.
    static void destroy(Mat& mat) noexcept;
};

/// @brief Throw std::runtime_error if a PETSc error code indicates failure.
void petsc_check(PetscErrorCode ierr, const char* action);

} // namespace gss
