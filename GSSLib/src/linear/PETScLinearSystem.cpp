#include "gss/linear/PETScLinearSystem.hpp"

#include <sstream>

namespace gss {

void petsc_check(PetscErrorCode ierr, const char* action) {
    if (ierr != 0) {
        std::ostringstream oss;
        oss << "PETSc failure while attempting to " << action << ": error code " << ierr;
        throw std::runtime_error(oss.str());
    }
}

Vec PETScLinearSystem::make_vec(Index global_size, MPI_Comm comm) {
    Vec vec{};
    petsc_check(VecCreateMPI(comm, PETSC_DECIDE, global_size, &vec), "create distributed vector");
    petsc_check(VecSet(vec, 0.0), "initialize vector");
    return vec;
}

Vec PETScLinearSystem::make_vec(std::span<const Real> values, MPI_Comm comm) {
    Vec vec = make_vec(static_cast<Index>(values.size()), comm);
    PetscInt rstart = 0;
    PetscInt rend = 0;
    petsc_check(VecGetOwnershipRange(vec, &rstart, &rend), "get vector ownership range");
    for (PetscInt i = rstart; i < rend; ++i) {
        const Real v = values[static_cast<std::size_t>(i)];
        petsc_check(VecSetValue(vec, i, v, INSERT_VALUES), "set vector value");
    }
    petsc_check(VecAssemblyBegin(vec), "begin vector assembly");
    petsc_check(VecAssemblyEnd(vec), "end vector assembly");
    return vec;
}

std::vector<Real> PETScLinearSystem::to_std_vector(Vec vec) {
    PetscInt n = 0;
    petsc_check(VecGetSize(vec, &n), "get vector size");
    std::vector<Real> out(static_cast<std::size_t>(n), 0.0);
    PetscInt rstart = 0;
    PetscInt rend = 0;
    petsc_check(VecGetOwnershipRange(vec, &rstart, &rend), "get vector ownership range");
    const PetscScalar* local = nullptr;
    petsc_check(VecGetArrayRead(vec, &local), "read vector array");
    for (PetscInt i = rstart; i < rend; ++i) {
        out[static_cast<std::size_t>(i)] = local[i - rstart];
    }
    petsc_check(VecRestoreArrayRead(vec, &local), "restore vector array");
    petsc_check(MPI_Allreduce(MPI_IN_PLACE, out.data(), n, MPIU_SCALAR, MPI_SUM, PetscObjectComm(reinterpret_cast<PetscObject>(vec))),
                "gather distributed vector");
    return out;
}

void PETScLinearSystem::destroy(Vec& vec) noexcept {
    if (vec != nullptr) {
        VecDestroy(&vec);
        vec = nullptr;
    }
}

void PETScLinearSystem::destroy(Mat& mat) noexcept {
    if (mat != nullptr) {
        MatDestroy(&mat);
        mat = nullptr;
    }
}

} // namespace gss
