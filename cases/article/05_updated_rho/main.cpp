#include "../support/official.hpp"
int main() {
  std::filesystem::current_path(ARTICLE_CASE_DIR);
  return p3official::run_program(p3official::Strategy::MixedRhoEveryIteration, "05_updated_rho");
}
