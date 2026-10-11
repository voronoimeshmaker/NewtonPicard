#ifdef GSMAKER_USE_SLEPC

#include "gsmaker/analysis/spectral/SlepcEigenSolverController.hpp"

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>


#include "gsmaker/analysis/spectral/SlepcContext.hpp"
#include "gsmaker/linear_system/PetscContext.hpp"
#include "gsmaker/linear_system/PetscMatrixUpdate.hpp"

namespace gss {
namespace {

class VecGuard final {
public:
    explicit VecGuard(Vec value = nullptr)
        : value_(value) {}

    ~VecGuard() noexcept {
        if (value_ != nullptr) {
            (void)VecDestroy(&value_);
        }
    }

    VecGuard(const VecGuard&) = delete;
    VecGuard& operator=(const VecGuard&) = delete;

    [[nodiscard]]
    Vec get() const noexcept {
        return value_;
    }

private:
    Vec value_{nullptr};
};


[[nodiscard]]
std::string current_eps_type(EPS eps) {

    EPSType type = nullptr;

    petsc_check(
        EPSGetType(eps, &type),
        "query active EPS type");

    return type != nullptr
        ? std::string{type}
        : std::string{};
}


[[nodiscard]]
std::vector<std::complex<PetscReal>> extract_eigenvector(
    Vec real_part,
    Vec imaginary_part) {

    PetscInt size = 0;

    petsc_check(
        VecGetSize(real_part, &size),
        "query SLEPc eigenvector size");

    std::vector<std::complex<PetscReal>> result(
        static_cast<std::size_t>(size));

    const PetscScalar* real_data = nullptr;

    petsc_check(
        VecGetArrayRead(real_part, &real_data),
        "access SLEPc eigenvector values");

#if defined(PETSC_USE_COMPLEX)

    try {
        for (PetscInt i = 0; i < size; ++i) {
            const PetscScalar value =
                real_data[static_cast<std::size_t>(i)];

            result[static_cast<std::size_t>(i)] = {
                PetscRealPart(value),
                PetscImaginaryPart(value)};
        }
    } catch (...) {
        (void)VecRestoreArrayRead(real_part, &real_data);
        throw;
    }

    petsc_check(
        VecRestoreArrayRead(real_part, &real_data),
        "restore SLEPc eigenvector values");

#else

    if (imaginary_part == nullptr) {
        (void)VecRestoreArrayRead(real_part, &real_data);
        throw std::invalid_argument(
            "real-scalar SLEPc extraction requires an imaginary eigenvector Vec");
    }

    const PetscScalar* imaginary_data = nullptr;

    petsc_check(
        VecGetArrayRead(imaginary_part, &imaginary_data),
        "access SLEPc imaginary eigenvector values");

    try {
        for (PetscInt i = 0; i < size; ++i) {
            const auto index = static_cast<std::size_t>(i);
            result[index] = {
                static_cast<PetscReal>(real_data[index]),
                static_cast<PetscReal>(imaginary_data[index])};
        }
    } catch (...) {
        (void)VecRestoreArrayRead(imaginary_part, &imaginary_data);
        (void)VecRestoreArrayRead(real_part, &real_data);
        throw;
    }

    petsc_check(
        VecRestoreArrayRead(imaginary_part, &imaginary_data),
        "restore SLEPc imaginary eigenvector values");

    petsc_check(
        VecRestoreArrayRead(real_part, &real_data),
        "restore SLEPc eigenvector values");

#endif

    return result;
}

} // namespace


SlepcEigenSolverController::
SlepcEigenSolverController(
    const CsrMatrix& matrix,
    SlepcEigenConfig config)
    : config_(std::move(config)) {

    if (!SlepcContext::initialized()) {
        throw std::runtime_error(
            "SLEPc must be initialized before constructing SlepcEigenSolverController");
    }

    config_.validate();
    matrix.validate();

    if (matrix.rows() <= 1) {
        throw std::invalid_argument(
            "SlepcEigenSolverController requires a matrix of order at least two");
    }

    if (config_.requested_eigenpairs >= matrix.rows()) {
        throw std::invalid_argument(
            "requested_eigenpairs must be smaller than the matrix order for the iterative EPS gate");
    }

    try {
        matrix_ = PetscMatrixUpdate::create(matrix);

        petsc_check(
            EPSCreate(PETSC_COMM_SELF, &eps_),
            "create SLEPc EPS context");

        configure_eps();

    } catch (...) {
        if (eps_ != nullptr) {
            (void)EPSDestroy(&eps_);
        }

        if (matrix_ != nullptr) {
            (void)MatDestroy(&matrix_);
        }

        throw;
    }
}


SlepcEigenSolverController::~SlepcEigenSolverController() noexcept {

    if (eps_ != nullptr) {
        (void)EPSDestroy(&eps_);
    }

    if (matrix_ != nullptr) {
        (void)MatDestroy(&matrix_);
    }
}


void SlepcEigenSolverController::configure_eps() {

    petsc_check(
        EPSSetOperators(eps_, matrix_, nullptr),
        "set EPS operators");

    petsc_check(
        EPSSetProblemType(eps_, config_.problem_type),
        "set EPS problem type");

    petsc_check(
        EPSSetType(eps_, config_.eps_type.c_str()),
        "set EPS solver type");

    petsc_check(
        EPSSetDimensions(
            eps_,
            config_.requested_eigenpairs,
            PETSC_DEFAULT,
            PETSC_DEFAULT),
        "set EPS dimensions");

    petsc_check(
        EPSSetTolerances(
            eps_,
            config_.tolerance,
            config_.max_iterations),
        "set EPS tolerances");

    petsc_check(
        EPSSetWhichEigenpairs(eps_, config_.which),
        "set EPS eigenpair selection");

    if (!config_.options_prefix.empty()) {
        petsc_check(
            EPSSetOptionsPrefix(
                eps_,
                config_.options_prefix.c_str()),
            "set EPS options prefix");
    }

    if (config_.set_from_options) {
        petsc_check(
            EPSSetFromOptions(eps_),
            "apply EPS options database overrides");
    }
}


SlepcEigenResult SlepcEigenSolverController::solve() {

    const auto begin = std::chrono::steady_clock::now();

    petsc_check(
        EPSSolve(eps_),
        "solve eigenvalue problem");

    const auto end = std::chrono::steady_clock::now();

    SlepcEigenResult result;

    result.requested_eigenpairs =
        config_.requested_eigenpairs;

    result.solve_seconds =
        std::chrono::duration<double>(end - begin).count();

    petsc_check(
        EPSGetConverged(eps_, &result.converged_eigenpairs),
        "query number of converged EPS eigenpairs");

    petsc_check(
        EPSGetIterationNumber(eps_, &result.iterations),
        "query EPS iteration count");

    petsc_check(
        EPSGetConvergedReason(eps_, &result.reason),
        "query EPS converged reason");

    result.eps_type = current_eps_type(eps_);

    const PetscInt extract_count = std::min(
        result.converged_eigenpairs,
        config_.requested_eigenpairs);

    Vec real_vector = nullptr;
    Vec imaginary_vector = nullptr;

    if (config_.compute_eigenvectors && extract_count > 0) {
        petsc_check(
            MatCreateVecs(matrix_, &real_vector, nullptr),
            "create EPS real eigenvector work vector");

#if !defined(PETSC_USE_COMPLEX)
        petsc_check(
            VecDuplicate(real_vector, &imaginary_vector),
            "create EPS imaginary eigenvector work vector");
#endif
    }

    VecGuard real_guard(real_vector);
    VecGuard imaginary_guard(imaginary_vector);

    result.eigenpairs.reserve(
        static_cast<std::size_t>(extract_count));

    for (PetscInt i = 0; i < extract_count; ++i) {
        PetscScalar eigen_real{};

#if defined(PETSC_USE_COMPLEX)
        petsc_check(
            EPSGetEigenpair(
                eps_,
                i,
                &eigen_real,
                nullptr,
                config_.compute_eigenvectors ? real_guard.get() : nullptr,
                nullptr),
            "extract EPS eigenpair");
#else
        PetscScalar eigen_imaginary{};

        petsc_check(
            EPSGetEigenpair(
                eps_,
                i,
                &eigen_real,
                &eigen_imaginary,
                config_.compute_eigenvectors ? real_guard.get() : nullptr,
                config_.compute_eigenvectors ? imaginary_guard.get() : nullptr),
            "extract EPS eigenpair");
#endif

        SlepcEigenpair eigenpair;

#if defined(PETSC_USE_COMPLEX)
        eigenpair.eigenvalue = {
            PetscRealPart(eigen_real),
            PetscImaginaryPart(eigen_real)};
#else
        eigenpair.eigenvalue = {
            static_cast<PetscReal>(eigen_real),
            static_cast<PetscReal>(eigen_imaginary)};
#endif

        if (config_.compute_eigenvectors) {
            eigenpair.eigenvector = extract_eigenvector(
                real_guard.get(),
                imaginary_guard.get());
        }

        petsc_check(
            EPSComputeError(
                eps_,
                i,
                EPS_ERROR_ABSOLUTE,
                &eigenpair.absolute_error),
            "compute EPS absolute eigenpair error");

        petsc_check(
            EPSComputeError(
                eps_,
                i,
                EPS_ERROR_RELATIVE,
                &eigenpair.relative_error),
            "compute EPS relative eigenpair error");

        petsc_check(
            EPSComputeError(
                eps_,
                i,
                EPS_ERROR_BACKWARD,
                &eigenpair.backward_error),
            "compute EPS backward eigenpair error");

        result.eigenpairs.push_back(
            std::move(eigenpair));
    }

    return result;
}

} // namespace gss

#endif // GSMAKER_USE_SLEPC
