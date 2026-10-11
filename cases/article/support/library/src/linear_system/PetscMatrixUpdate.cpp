#include "gsmaker/linear_system/PetscMatrixUpdate.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "gsmaker/linear_system/PetscContext.hpp"


namespace gss {
namespace {

void insert_matrix_values(
    Mat target,
    const CsrMatrix& matrix) {

    const auto row_start =
        matrix.row_start();

    const auto column =
        matrix.column();

    const auto value =
        matrix.value();

    for (Index row = 0;
         row < matrix.rows();
         ++row) {

        const auto rs =
            static_cast<std::size_t>(
                row);

        const Index begin =
            row_start[rs];

        const Index end =
            row_start[rs + 1U];

        const PetscInt count =
            static_cast<PetscInt>(
                end - begin);

        if (count == 0) {
            continue;
        }

        const PetscInt petsc_row =
            static_cast<PetscInt>(
                row);

        petsc_check(
            MatSetValues(
                target,
                1,
                &petsc_row,
                count,
                column.data() +
                    static_cast<std::size_t>(
                        begin),
                value.data() +
                    static_cast<std::size_t>(
                        begin),
                INSERT_VALUES),
            "insert CSR row into PETSc matrix");
    }

    petsc_check(
        MatAssemblyBegin(
            target,
            MAT_FINAL_ASSEMBLY),
        "begin PETSc matrix assembly");

    petsc_check(
        MatAssemblyEnd(
            target,
            MAT_FINAL_ASSEMBLY),
        "finish PETSc matrix assembly");
}

} // namespace


Mat PetscMatrixUpdate::create(
    const CsrMatrix& matrix) {

    if (!PetscContext::initialized()) {
        throw std::runtime_error(
            "PETSc must be initialized before creating a PETSc matrix");
    }

    matrix.validate();

    const PetscInt n =
        static_cast<PetscInt>(
            matrix.rows());

    std::vector<PetscInt>
        nonzeros_per_row(
            static_cast<std::size_t>(
                n),
            PetscInt{0});

    const auto row_start =
        matrix.row_start();

    for (PetscInt row = 0;
         row < n;
         ++row) {

        const auto rs =
            static_cast<std::size_t>(
                row);

        nonzeros_per_row[rs] =
            static_cast<PetscInt>(
                row_start[rs + 1U] -
                row_start[rs]);
    }

    Mat result = nullptr;

    petsc_check(
        MatCreateSeqAIJ(
            PETSC_COMM_SELF,
            n,
            n,
            0,
            nonzeros_per_row.data(),
            &result),
        "create sequential AIJ matrix");

    try {

        insert_matrix_values(
            result,
            matrix);

        petsc_check(
            MatSetOption(
                result,
                MAT_NEW_NONZERO_LOCATION_ERR,
                PETSC_TRUE),
            "lock PETSc matrix sparsity pattern");

    } catch (...) {

        if (result != nullptr) {
            (void)MatDestroy(
                &result);
        }

        throw;
    }

    return result;
}


void PetscMatrixUpdate::update(
    Mat matrix,
    const CsrMatrix& values) {

    if (matrix == nullptr) {
        throw std::invalid_argument(
            "PetscMatrixUpdate received a null PETSc matrix");
    }

    values.validate();

    PetscInt rows = 0;
    PetscInt columns = 0;

    petsc_check(
        MatGetSize(
            matrix,
            &rows,
            &columns),
        "query PETSc matrix dimensions");

    if (rows !=
            static_cast<PetscInt>(
                values.rows()) ||
        columns != rows) {

        throw std::invalid_argument(
            "PetscMatrixUpdate matrix dimension mismatch");
    }

    petsc_check(
        MatZeroEntries(
            matrix),
        "zero PETSc matrix entries");

    insert_matrix_values(
        matrix,
        values);
}

} // namespace gss
