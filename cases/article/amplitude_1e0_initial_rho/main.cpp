#include "../support/complementary.hpp"
int main() {
  std::filesystem::current_path(ARTICLE_CASE_DIR);
  return p3official::run_program(p3official::Strategy::MixedRhoInitial, "amplitude_1e0_initial_rho");
}
