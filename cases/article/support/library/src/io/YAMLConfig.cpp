#include "gsmaker/io/YAMLConfig.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace gss {

namespace {

[[nodiscard]]
YAML::Node section(
    const YAML::Node& root,
    const char* key) {

    if (!root || !root.IsMap()) {
        return {};
    }

    const YAML::Node node = root[key];

    if (!node || !node.IsMap()) {
        return {};
    }

    return node;
}


template<typename T>
[[nodiscard]]
T require_as(
    const YAML::Node& parent,
    const char* key,
    const char* path) {

    try {
        if (!parent || !parent.IsMap()) {
            throw std::runtime_error("");
        }

        const YAML::Node value = parent[key];

        if (!value) {
            throw std::runtime_error("");
        }

        return value.as<T>();

    } catch (const std::exception&) {
        throw std::runtime_error(
            std::string("invalid or missing YAML field: ") + path);
    }
}


template<typename T>
[[nodiscard]]
T optional_as(
    const YAML::Node& parent,
    const char* key,
    const T& fallback,
    const char* path) {

    if (!parent || !parent.IsMap()) {
        return fallback;
    }

    const YAML::Node value = parent[key];

    if (!value) {
        return fallback;
    }

    try {
        return value.as<T>();
    } catch (const std::exception&) {
        throw std::runtime_error(
            std::string("invalid YAML field: ") + path);
    }
}

} // namespace


Config YAMLConfig::load(
    const std::filesystem::path& path) {

    const YAML::Node root =
        YAML::LoadFile(path.string());

    const YAML::Node mesh =
        section(root, "mesh");

    const YAML::Node solver =
        section(root, "solver");

    const YAML::Node newton =
        section(root, "newton");

    const YAML::Node petsc =
        section(root, "petsc");

    Config config{};


    config.mesh.type =
        require_as<std::string>(
            mesh,
            "type",
            "mesh.type");

    config.mesh.n_generators =
        require_as<Index>(
            mesh,
            "n_generators",
            "mesh.n_generators");

    config.mesh.n_per_direction =
        optional_as(
            mesh,
            "n_per_direction",
            config.mesh.n_per_direction,
            "mesh.n_per_direction");

    config.mesh.n_per_direction_y =
        optional_as(
            mesh,
            "n_per_direction_y",
            config.mesh.n_per_direction_y,
            "mesh.n_per_direction_y");

    config.mesh.random_seed =
        optional_as(
            mesh,
            "random_seed",
            config.mesh.random_seed,
            "mesh.random_seed");

    const auto domain =
        require_as<std::vector<Real>>(
            mesh,
            "domain",
            "mesh.domain");

    if (domain.size() != 4) {
        throw std::runtime_error(
            "invalid YAML field: mesh.domain must contain four values");
    }

    std::copy(
        domain.begin(),
        domain.end(),
        config.mesh.domain.begin());

    config.mesh.hilbert_reorder =
        optional_as(
            mesh,
            "hilbert_reorder",
            config.mesh.hilbert_reorder,
            "mesh.hilbert_reorder");

    config.mesh.reordering =
        optional_as(
            mesh,
            "reordering",
            config.mesh.reordering,
            "mesh.reordering");


    config.solver.nonlinear =
        require_as<std::string>(
            solver,
            "nonlinear",
            "solver.nonlinear");

    config.solver.linear =
        require_as<std::string>(
            solver,
            "linear",
            "solver.linear");

    config.solver.max_iter =
        require_as<Index>(
            solver,
            "max_iter",
            "solver.max_iter");

    config.solver.tolerance_relative =
        require_as<Real>(
            solver,
            "tolerance_relative",
            "solver.tolerance_relative");

    config.solver.tolerance_absolute =
        require_as<Real>(
            solver,
            "tolerance_absolute",
            "solver.tolerance_absolute");

    config.solver.relaxation_omega =
        optional_as(
            solver,
            "relaxation_omega",
            config.solver.relaxation_omega,
            "solver.relaxation_omega");

    config.solver.rhs_constant =
        optional_as(
            solver,
            "rhs_constant",
            config.solver.rhs_constant,
            "solver.rhs_constant");

    config.solver.source_regularisation_delta =
        optional_as(
            solver,
            "source_regularisation_delta",
            config.solver.source_regularisation_delta,
            "solver.source_regularisation_delta");

    config.solver.diagonal_shift =
        optional_as(
            solver,
            "diagonal_shift",
            config.solver.diagonal_shift,
            "solver.diagonal_shift");

    config.solver.initial_guess =
        optional_as(
            solver,
            "initial_guess",
            config.solver.initial_guess,
            "solver.initial_guess");

    config.solver.perturbation_amplitude =
        optional_as(
            solver,
            "perturbation_amplitude",
            config.solver.perturbation_amplitude,
            "solver.perturbation_amplitude");

    config.solver.perturbation_type =
        optional_as(
            solver,
            "perturbation_type",
            config.solver.perturbation_type,
            "solver.perturbation_type");


    config.newton.jacobian =
        optional_as(
            newton,
            "jacobian",
            config.newton.jacobian,
            "newton.jacobian");

    config.newton.fd_epsilon =
        optional_as(
            newton,
            "fd_epsilon",
            config.newton.fd_epsilon,
            "newton.fd_epsilon");

    config.newton.linesearch =
        optional_as(
            newton,
            "linesearch",
            config.newton.linesearch,
            "newton.linesearch");

    config.newton.linesearch_max_iter =
        optional_as(
            newton,
            "linesearch_max_iter",
            config.newton.linesearch_max_iter,
            "newton.linesearch_max_iter");

    config.newton.linesearch_reduction =
        optional_as(
            newton,
            "linesearch_reduction",
            config.newton.linesearch_reduction,
            "newton.linesearch_reduction");


    config.petsc.ksp_type =
        optional_as(
            petsc,
            "ksp_type",
            config.petsc.ksp_type,
            "petsc.ksp_type");

    config.petsc.pc_type =
        optional_as(
            petsc,
            "pc_type",
            config.petsc.pc_type,
            "petsc.pc_type");

    config.petsc.mat_solver_type =
        optional_as(
            petsc,
            "mat_solver_type",
            config.petsc.mat_solver_type,
            "petsc.mat_solver_type");

    config.petsc.mumps_reorder =
        optional_as(
            petsc,
            "mumps_reorder",
            config.petsc.mumps_reorder,
            "petsc.mumps_reorder");

    config.petsc.mumps_icntl_14 =
        optional_as(
            petsc,
            "mumps_icntl_14",
            config.petsc.mumps_icntl_14,
            "petsc.mumps_icntl_14");

    return config;
}


void YAMLConfig::write(
    const Config& config,
    const std::filesystem::path& path) {

    YAML::Emitter out;

    out << YAML::BeginMap;


    out
        << YAML::Key << "mesh"
        << YAML::Value << YAML::BeginMap

        << YAML::Key << "type"
        << YAML::Value << config.mesh.type

        << YAML::Key << "n_generators"
        << YAML::Value << config.mesh.n_generators

        << YAML::Key << "n_per_direction"
        << YAML::Value << config.mesh.n_per_direction

        << YAML::Key << "n_per_direction_y"
        << YAML::Value << config.mesh.n_per_direction_y

        << YAML::Key << "random_seed"
        << YAML::Value << config.mesh.random_seed

        << YAML::Key << "domain"
        << YAML::Value
        << YAML::Flow
        << std::vector<Real>{
            config.mesh.domain.begin(),
            config.mesh.domain.end()
        }

        << YAML::Key << "hilbert_reorder"
        << YAML::Value << config.mesh.hilbert_reorder

        << YAML::Key << "reordering"
        << YAML::Value << config.mesh.reordering

        << YAML::EndMap;


    out
        << YAML::Key << "solver"
        << YAML::Value << YAML::BeginMap

        << YAML::Key << "nonlinear"
        << YAML::Value << config.solver.nonlinear

        << YAML::Key << "linear"
        << YAML::Value << config.solver.linear

        << YAML::Key << "max_iter"
        << YAML::Value << config.solver.max_iter

        << YAML::Key << "tolerance_relative"
        << YAML::Value << config.solver.tolerance_relative

        << YAML::Key << "tolerance_absolute"
        << YAML::Value << config.solver.tolerance_absolute

        << YAML::Key << "relaxation_omega"
        << YAML::Value << config.solver.relaxation_omega

        << YAML::Key << "rhs_constant"
        << YAML::Value << config.solver.rhs_constant

        << YAML::Key << "source_regularisation_delta"
        << YAML::Value
        << config.solver.source_regularisation_delta

        << YAML::Key << "diagonal_shift"
        << YAML::Value << config.solver.diagonal_shift

        << YAML::Key << "initial_guess"
        << YAML::Value << config.solver.initial_guess

        << YAML::Key << "perturbation_amplitude"
        << YAML::Value
        << config.solver.perturbation_amplitude

        << YAML::Key << "perturbation_type"
        << YAML::Value
        << config.solver.perturbation_type

        << YAML::EndMap;


    out
        << YAML::Key << "newton"
        << YAML::Value << YAML::BeginMap

        << YAML::Key << "jacobian"
        << YAML::Value << config.newton.jacobian

        << YAML::Key << "fd_epsilon"
        << YAML::Value << config.newton.fd_epsilon

        << YAML::Key << "linesearch"
        << YAML::Value << config.newton.linesearch

        << YAML::Key << "linesearch_max_iter"
        << YAML::Value
        << config.newton.linesearch_max_iter

        << YAML::Key << "linesearch_reduction"
        << YAML::Value
        << config.newton.linesearch_reduction

        << YAML::EndMap;


    out
        << YAML::Key << "petsc"
        << YAML::Value << YAML::BeginMap

        << YAML::Key << "ksp_type"
        << YAML::Value << config.petsc.ksp_type

        << YAML::Key << "pc_type"
        << YAML::Value << config.petsc.pc_type

        << YAML::Key << "mat_solver_type"
        << YAML::Value
        << config.petsc.mat_solver_type

        << YAML::Key << "mumps_reorder"
        << YAML::Value << config.petsc.mumps_reorder

        << YAML::Key << "mumps_icntl_14"
        << YAML::Value << config.petsc.mumps_icntl_14

        << YAML::EndMap;


    out << YAML::EndMap;


    std::ofstream file(path);

    if (!file) {
        throw std::runtime_error(
            "failed to open YAML output: " +
            path.string());
    }

    file << out.c_str() << '\n';
}

} // namespace gss