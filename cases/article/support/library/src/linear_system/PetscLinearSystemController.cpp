#include "gsmaker/linear_system/PetscLinearSystemController.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "gsmaker/linear_system/PetscContext.hpp"
#include "gsmaker/linear_system/PetscMatrixUpdate.hpp"




namespace gss {
namespace {

class VecGuard final {
public:
    VecGuard() = default;

    explicit VecGuard(
        Vec value)
        : value_(value) {}

    ~VecGuard() noexcept {
        if (value_ != nullptr) {
            (void)VecDestroy(
                &value_);
        }
    }

    VecGuard(
        const VecGuard&) = delete;

    VecGuard& operator=(
        const VecGuard&) = delete;

    [[nodiscard]]
    Vec get() const noexcept {
        return value_;
    }

private:
    Vec value_{nullptr};
};


[[nodiscard]]
Vec make_vector(
    Index size) {

    Vec vector = nullptr;

    petsc_check(
        VecCreateSeq(
            PETSC_COMM_SELF,
            static_cast<PetscInt>(
                size),
            &vector),
        "create sequential PETSc vector");

    return vector;
}


void assign_vector(
    Vec target,
    std::span<const Real> values) {

    PetscInt size = 0;

    petsc_check(
        VecGetSize(
            target,
            &size),
        "query PETSc vector size");

    if (static_cast<std::size_t>(
            size) !=
        values.size()) {

        throw std::invalid_argument(
            "PETSc vector assignment size mismatch");
    }

    PetscScalar* data = nullptr;

    petsc_check(
        VecGetArray(
            target,
            &data),
        "access writable PETSc vector array");

    try {

        std::copy(
            values.begin(),
            values.end(),
            data);

    } catch (...) {

        (void)VecRestoreArray(
            target,
            &data);

        throw;
    }

    petsc_check(
        VecRestoreArray(
            target,
            &data),
        "restore writable PETSc vector array");
}


[[nodiscard]]
std::vector<Real> extract_vector(
    Vec source) {

    PetscInt size = 0;

    petsc_check(
        VecGetSize(
            source,
            &size),
        "query PETSc solution size");

    std::vector<Real> result(
        static_cast<std::size_t>(
            size));

    const PetscScalar* data =
        nullptr;

    petsc_check(
        VecGetArrayRead(
            source,
            &data),
        "access PETSc solution array");

    try {

        std::copy(
            data,
            data +
                static_cast<std::size_t>(
                    size),
            result.begin());

    } catch (...) {

        (void)VecRestoreArrayRead(
            source,
            &data);

        throw;
    }

    petsc_check(
        VecRestoreArrayRead(
            source,
            &data),
        "restore PETSc solution array");

    return result;
}


[[nodiscard]]
std::string current_ksp_type(
    KSP ksp) {

    KSPType type = nullptr;

    petsc_check(
        KSPGetType(
            ksp,
            &type),
        "query active KSP type");

    return type != nullptr
        ? std::string{type}
        : std::string{};
}




void apply_known_matrix_properties(
    Mat target,
    const CsrMatrix& source,
    const PetscLinearConfig& config) {

    const auto& properties = source.properties();

    const bool symmetric =
        properties.symmetric ||
        config.matrix_symmetric_hint;

    const bool structurally_symmetric =
        properties.structurally_symmetric ||
        symmetric;

    const bool spd =
        properties.positive_definite ||
        config.matrix_spd_hint;

    if (structurally_symmetric) {
        petsc_check(
            MatSetOption(
                target,
                MAT_STRUCTURALLY_SYMMETRIC,
                PETSC_TRUE),
            "declare PETSc matrix structurally symmetric");
    }

    if (symmetric) {
        petsc_check(
            MatSetOption(
                target,
                MAT_SYMMETRIC,
                PETSC_TRUE),
            "declare PETSc matrix symmetric");
    }

    if (spd) {
        petsc_check(
            MatSetOption(
                target,
                MAT_SPD,
                PETSC_TRUE),
            "declare PETSc matrix SPD");
    }
}

[[nodiscard]]
std::string current_pc_type(
    KSP ksp) {

    PC pc = nullptr;

    petsc_check(
        KSPGetPC(
            ksp,
            &pc),
        "query KSP preconditioner");

    PCType type = nullptr;

    petsc_check(
        PCGetType(
            pc,
            &type),
        "query active PC type");

    return type != nullptr
        ? std::string{type}
        : std::string{};
}

} // namespace


PetscLinearSystemController::
PetscLinearSystemController(
    const CsrMatrix& matrix,
    PetscLinearConfig config)
    : config_(std::move(config)),
      rows_(matrix.rows()) {

    if (!PetscContext::initialized()) {
        throw std::runtime_error(
            "PETSc must be initialized before constructing PetscLinearSystemController");
    }

    if (rows_ <= 0) {
        throw std::invalid_argument(
            "PetscLinearSystemController requires a non-empty matrix");
    }

    config_.validate();

    const auto begin =
        std::chrono::steady_clock::now();

    try {

        matrix_ =
            PetscMatrixUpdate::create(
                matrix);

        // Prefer mathematical properties carried by the matrix itself.
        // Legacy caller hints remain supported for matrices assembled outside
        // GSMakerLib or for properties (such as SPD) that depend on the
        // particular boundary-value problem.
        apply_known_matrix_properties(
            matrix_,
            matrix,
            config_);

        petsc_check(
            KSPCreate(
                PETSC_COMM_SELF,
                &ksp_),
            "create PETSc KSP");

        petsc_check(
            KSPSetOperators(
                ksp_,
                matrix_,
                matrix_),
            "attach operator to PETSc KSP");

        configure_ksp();

        petsc_check(
            KSPSetUp(
                ksp_),
            "set up PETSc KSP");

    } catch (...) {

        if (ksp_ != nullptr) {
            (void)KSPDestroy(
                &ksp_);
        }

        if (matrix_ != nullptr) {
            (void)MatDestroy(
                &matrix_);
        }

        throw;
    }

    setup_seconds_ =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            begin)
            .count();
}


PetscLinearSystemController::
~PetscLinearSystemController() noexcept {

    if (ksp_ != nullptr) {
        (void)KSPDestroy(
            &ksp_);
    }

    if (matrix_ != nullptr) {
        (void)MatDestroy(
            &matrix_);
    }
}


void PetscLinearSystemController::
configure_ksp() {

    petsc_check(
        KSPSetType(
            ksp_,
            config_.ksp_type.c_str()),
        "set PETSc KSP type");

    PC pc = nullptr;

    petsc_check(
        KSPGetPC(
            ksp_,
            &pc),
        "query KSP preconditioner");

    petsc_check(
        PCSetType(
            pc,
            config_.pc_type.c_str()),
        "set PETSc PC type");

    if (!config_.mat_solver_type.empty() &&
        (config_.pc_type == "lu" ||
         config_.pc_type == "cholesky")) {

        petsc_check(
            PCFactorSetMatSolverType(
                pc,
                config_.mat_solver_type.c_str()),
            "set PETSc factorization package");
    }

    petsc_check(
        KSPSetTolerances(
            ksp_,
            config_.relative_tolerance,
            config_.absolute_tolerance,
            config_.divergence_tolerance,
            config_.max_iterations),
        "set PETSc KSP tolerances");

    if (!config_.options_prefix.empty()) {

        petsc_check(
            KSPSetOptionsPrefix(
                ksp_,
                config_.options_prefix.c_str()),
            "set PETSc KSP options prefix");
    }

    if (config_.set_from_options) {

        // This comes after the programmatic defaults deliberately: PETSc's
        // own options database remains the final authority when the user
        // supplies an override.
        petsc_check(
            KSPSetFromOptions(
                ksp_),
            "apply PETSc options database to KSP");
    }
}


PetscLinearResult
PetscLinearSystemController::solve(
    std::span<const Real> rhs,
    std::span<const Real> initial_guess) {

    if (rhs.size() !=
        static_cast<std::size_t>(
            rows_)) {

        throw std::invalid_argument(
            "PetscLinearSystemController RHS size mismatch");
    }

    if (!initial_guess.empty() &&
        initial_guess.size() !=
            rhs.size()) {

        throw std::invalid_argument(
            "PetscLinearSystemController initial guess size mismatch");
    }

    for (const auto value :
         rhs) {

        if (!std::isfinite(
                static_cast<double>(
                    value))) {

            throw std::invalid_argument(
                "PetscLinearSystemController RHS contains a non-finite value");
        }
    }

    VecGuard b{
        make_vector(
            rows_)};

    VecGuard x{
        make_vector(
            rows_)};

    assign_vector(
        b.get(),
        rhs);

    if (initial_guess.empty()) {

        petsc_check(
            VecSet(
                x.get(),
                PetscScalar{0}),
            "zero PETSc initial guess");

        petsc_check(
            KSPSetInitialGuessNonzero(
                ksp_,
                PETSC_FALSE),
            "declare zero PETSc initial guess");

    } else {

        assign_vector(
            x.get(),
            initial_guess);

        petsc_check(
            KSPSetInitialGuessNonzero(
                ksp_,
                PETSC_TRUE),
            "declare non-zero PETSc initial guess");
    }

    const auto solve_begin =
        std::chrono::steady_clock::now();

    petsc_check(
        KSPSolve(
            ksp_,
            b.get(),
            x.get()),
        "solve PETSc linear system");

    const double solve_seconds =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            solve_begin)
            .count();

    PetscInt iterations = 0;
    PetscReal ksp_residual = 0;

    KSPConvergedReason reason =
        KSP_CONVERGED_ITERATING;

    petsc_check(
        KSPGetIterationNumber(
            ksp_,
            &iterations),
        "query PETSc KSP iteration count");

    petsc_check(
        KSPGetResidualNorm(
            ksp_,
            &ksp_residual),
        "query PETSc KSP residual norm");

    petsc_check(
        KSPGetConvergedReason(
            ksp_,
            &reason),
        "query PETSc KSP convergence reason");

    const char* reason_text =
        nullptr;

    petsc_check(
        KSPGetConvergedReasonString(
            ksp_,
            &reason_text),
        "query PETSc KSP convergence reason string");

    Vec residual = nullptr;

    petsc_check(
        VecDuplicate(
            b.get(),
            &residual),
        "create PETSc true-residual vector");

    VecGuard residual_guard{
        residual};

    petsc_check(
        MatMult(
            matrix_,
            x.get(),
            residual_guard.get()),
        "evaluate PETSc true residual");

    petsc_check(
        VecAXPY(
            residual_guard.get(),
            PetscScalar{-1},
            b.get()),
        "subtract RHS from PETSc true residual");

    PetscReal true_residual = 0;
    PetscReal rhs_norm = 0;

    petsc_check(
        VecNorm(
            residual_guard.get(),
            NORM_2,
            &true_residual),
        "compute PETSc true residual norm");

    petsc_check(
        VecNorm(
            b.get(),
            NORM_2,
            &rhs_norm),
        "compute PETSc RHS norm");

    const PetscReal scale =
        std::max(
            rhs_norm,
            std::numeric_limits<PetscReal>::
                min());

    return PetscLinearResult{
        extract_vector(
            x.get()),
        static_cast<Index>(
            iterations),
        ksp_residual,
        true_residual,
        rhs_norm,
        true_residual / scale,
        reason,
        reason_text != nullptr
            ? std::string{reason_text}
            : std::string{},
        current_ksp_type(
            ksp_),
        current_pc_type(
            ksp_),
        setup_seconds_,
        solve_seconds};
}


void PetscLinearSystemController::
update_matrix(
    const CsrMatrix& matrix) {

    if (matrix.rows() != rows_) {
        throw std::invalid_argument(
            "PetscLinearSystemController matrix dimension changed");
    }

    const auto begin =
        std::chrono::steady_clock::now();

    PetscMatrixUpdate::update(
        matrix_,
        matrix);

    // PetscMatrixUpdate changes matrix values.  Reassert the properties of
    // the replacement CsrMatrix after assembly instead of assuming that the
    // previous operator's mathematical properties still apply.
    apply_known_matrix_properties(
        matrix_,
        matrix,
        config_);

    petsc_check(
        KSPSetOperators(
            ksp_,
            matrix_,
            matrix_),
        "reattach updated operator to PETSc KSP");

    petsc_check(
        KSPSetUp(
            ksp_),
        "set up PETSc KSP after matrix update");

    setup_seconds_ =
        std::chrono::duration<double>(
            std::chrono::steady_clock::now() -
            begin)
            .count();
}

} // namespace gss
