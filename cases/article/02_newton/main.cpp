#include "../support/official.hpp"
int main() {
  std::filesystem::current_path(ARTICLE_CASE_DIR);
  return p3official::run_program(p3official::Strategy::NewtonOnly, "02_newton");
}
