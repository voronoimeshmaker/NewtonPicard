#include "../support/spatial.cpp"
int main() {
  std::filesystem::current_path(ARTICLE_CASE_DIR);
  return spatial_main();
}
