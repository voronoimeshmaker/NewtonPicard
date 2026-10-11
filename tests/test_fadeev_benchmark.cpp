#include <gtest/gtest.h>

#include "gss/benchmark/FadeevSolution.hpp"

TEST(FadeevBenchmark, DerivativesAreFinite) {
    const gss::FadeevSolution f{};
    EXPECT_TRUE(std::isfinite(f.psi_exact(0.1, 0.2)));
    EXPECT_TRUE(std::isfinite(f.laplacian(0.1, 0.2)));
    EXPECT_FALSE(f.xpoint_locations().empty());
    EXPECT_FALSE(f.opoint_locations().empty());
}
