#include <gtest/gtest.h>

#include "gss/solvers/NewtonSolver.hpp"
#include "gss/linear/CholeskyBackend.hpp"

TEST(NewtonSolver, FactoryTypeIsConstructible) {
    gss::Config config{};
    EXPECT_NO_THROW((gss::NewtonSolver<gss::CholeskyBackend>{config}));
}
