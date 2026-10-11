#include "../support/complementary.hpp"
int main() {
  std::filesystem::current_path(ARTICLE_CASE_DIR);
  return p3official::run_program(p3official::Strategy::MixedRhoEveryIteration, "amplitude_1e-2_updated_rho");
}
