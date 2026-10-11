#include <gtest/gtest.h>

#include "gss/solvers/PicardSolver.hpp"
#include "gss/linear/CholeskyBackend.hpp"

TEST(PicardSolver, FactoryTypeIsConstructible) {
    gss::Config config{};
    EXPECT_NO_THROW((gss::PicardSolver<gss::CholeskyBackend>{config}));
}
