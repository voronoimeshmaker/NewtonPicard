#include "gss/linear/CholeskyBackend.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

[[nodiscard]] int mumps_ordering_code(std::string_view name) {
    constexpr std::array<std::pair<std::string_view, int>, 3> orderings{{
        {"amd", 0}, {"pord", 4}, {"metis", 5}}};
    const auto found = std::find_if(
        orderings.begin(), orderings.end(),
        [name](const auto& entry) { return entry.first == name; });
    if (found == orderings.end()) {
        throw std::invalid_argument(
            "Unknown MUMPS ordering; expected amd, pord or metis");
    }
    return found->second;
}

class ScopedPetscOption {
public:
    ScopedPetscOption(std::string key, std::string value)
        : key_(std::move(key)) {
        std::array<char, 64> previous{};
        PetscBool present = PETSC_FALSE;
        gss::petsc_check(
            PetscOptionsGetString(nullptr, nullptr, key_.c_str(), previous.data(),
                                  previous.size(), &present),
            "read existing PETSc option");
        if (present == PETSC_TRUE) {
            previous_ = previous.data();
        }
        gss::petsc_check(
            PetscOptionsSetValue(nullptr, key_.c_str(), value.c_str()),
            "set scoped PETSc option");
    }

    ~ScopedPetscOption() {
        if (previous_.has_value()) {
            PetscOptionsSetValue(nullptr, key_.c_str(), previous_->c_str());
        } else {
            PetscOptionsClearValue(nullptr, key_.c_str());
        }
    }

    ScopedPetscOption(const ScopedPetscOption&) = delete;
    ScopedPetscOption& operator=(const ScopedPetscOption&) = delete;

private:
    std::string key_;
    std::optional<std::string> previous_;
};

} // namespace

namespace gss {

CholeskyBackend::~CholeskyBackend() {
    if (ksp_ != nullptr) {
        KSPDestroy(&ksp_);
    }
}

CholeskyBackend::CholeskyBackend(CholeskyBackend&& other) noexcept
    : config_(std::move(other.config_)), ksp_(other.ksp_) {
    other.ksp_ = nullptr;
}

CholeskyBackend& CholeskyBackend::operator=(CholeskyBackend&& other) noexcept {
    if (this != &other) {
        if (ksp_ != nullptr) {
            KSPDestroy(&ksp_);
        }
        config_ = std::move(other.config_);
        ksp_ = other.ksp_;
        other.ksp_ = nullptr;
    }
    return *this;
}

void CholeskyBackend::factorize(Mat A) {
    if (ksp_ != nullptr) {
        KSPDestroy(&ksp_);
        ksp_ = nullptr;
    }
    MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(A));
    petsc_check(KSPCreate(comm, &ksp_), "create KSP");
    petsc_check(KSPSetOperators(ksp_, A, A), "set KSP operators");
    std::string ksp_type = config_.petsc.ksp_type;
    std::string pc_type = config_.petsc.pc_type;
    if (config_.solver.linear == "cholesky_sparse") {
        ksp_type = "preonly";
        pc_type = "cholesky";
    } else if (config_.solver.linear == "cg_jacobi") {
        ksp_type = "cg";
        pc_type = "jacobi";
    } else if (config_.solver.linear == "cg_gamg") {
        ksp_type = "cg";
        pc_type = "gamg";
    }
    petsc_check(KSPSetType(ksp_, ksp_type.c_str()), "set KSP type");
    PC pc{};
    petsc_check(KSPGetPC(ksp_, &pc), "get PC");
    petsc_check(PCSetType(pc, pc_type.c_str()), "set PC type");
    if (pc_type == "cholesky" || pc_type == "lu") {
        petsc_check(PCFactorSetMatSolverType(pc, config_.petsc.mat_solver_type.c_str()), "set matrix solver type");
    }
    petsc_check(KSPSetTolerances(
        ksp_,
        config_.solver.tolerance_relative * Real{0.1},
        config_.solver.tolerance_absolute * Real{0.1},
        PETSC_DEFAULT,
        PETSC_DEFAULT), "set KSP tolerances");
    std::optional<ScopedPetscOption> mumps_ordering;
    if ((pc_type == "cholesky" || pc_type == "lu")
        && config_.petsc.mat_solver_type == "mumps") {
        mumps_ordering.emplace(
            "-mat_mumps_icntl_7",
            std::to_string(mumps_ordering_code(config_.petsc.mumps_reorder)));
    }
    petsc_check(KSPSetFromOptions(ksp_), "apply KSP options");
    petsc_check(KSPSetUp(ksp_), "factorize matrix");
}

void CholeskyBackend::solve(Vec b, Vec x) {
    petsc_check(KSPSolve(ksp_, b, x), "solve linear system");
    PetscInt its = 0;
    petsc_check(KSPGetIterationNumber(ksp_, &its), "read KSP iteration count");
    last_iterations_ = its;
    KSPConvergedReason reason{};
    petsc_check(KSPGetConvergedReason(ksp_, &reason), "read KSP convergence reason");
    if (reason < 0) {
        throw std::runtime_error("PETSc KSP diverged");
    }
}

} // namespace gss
