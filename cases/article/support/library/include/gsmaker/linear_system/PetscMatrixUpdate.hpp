#pragma once

#include <petscmat.h>

#include "gsmaker/assembly/CsrMatrix.hpp"


namespace gss {

/// @brief Sequential PETSc MATSEQAIJ bridge for GSMakerLib CSR matrices.
///
/// This phase is intentionally serial and uses PETSC_COMM_SELF.  Parallel
/// distribution belongs to a later solver/backend phase.
class PetscMatrixUpdate final {
public:
    [[nodiscard]]
    static Mat create(
        const CsrMatrix& matrix);

    static void update(
        Mat matrix,
        const CsrMatrix& values);
};

} // namespace gss
