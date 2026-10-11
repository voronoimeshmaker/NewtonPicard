// SPDX-License-Identifier: GPL-3.0-or-later
/** @file test_config.cpp
 * @brief Configuration regressions independent of the build-directory depth.
 * @ingroup tests
 */
#include <gtest/gtest.h>

#include "gss/io/YAMLConfig.hpp"

namespace {

[[nodiscard]] std::filesystem::path default_config_path() {
    return std::filesystem::path(GSS_TEST_SOURCE_DIR).parent_path()/"config/default.yaml";
}

} // namespace

TEST(Config, DefaultsAreUsable) {
    const auto config = gss::YAMLConfig::load(default_config_path());
    EXPECT_EQ(config.solver.nonlinear, "picard");
    EXPECT_GT(config.mesh.n_generators, 0);
    EXPECT_GT(config.mesh.n_per_direction, 0);
    EXPECT_GT(config.solver.tolerance_relative, 0.0);
    EXPECT_DOUBLE_EQ(config.solver.source_regularisation_delta, 1.0e-10);
    EXPECT_EQ(config.petsc.mumps_reorder, "pord");
    EXPECT_DOUBLE_EQ(config.benchmark.fadeev_alpha, 0.225);
    EXPECT_DOUBLE_EQ(config.benchmark.fadeev_theta_degrees, 5.7);
}

TEST(Config, WriteReadRoundTrip) {
    gss::Config config{};
    const std::filesystem::path path = "config_roundtrip.yaml";
    gss::YAMLConfig::write(config, path);
    const auto parsed = gss::YAMLConfig::load(path);
    EXPECT_EQ(parsed.mesh.type, config.mesh.type);
    EXPECT_EQ(parsed.solver.linear, config.solver.linear);
    EXPECT_EQ(parsed.petsc.mumps_reorder, config.petsc.mumps_reorder);
    EXPECT_DOUBLE_EQ(parsed.solver.source_regularisation_delta, config.solver.source_regularisation_delta);
    EXPECT_DOUBLE_EQ(parsed.benchmark.fadeev_alpha, config.benchmark.fadeev_alpha);
    EXPECT_DOUBLE_EQ(parsed.benchmark.fadeev_theta_degrees, config.benchmark.fadeev_theta_degrees);
    std::filesystem::remove(path);
}
