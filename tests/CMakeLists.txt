file(GLOB TEST_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp")

foreach(test_source IN LISTS TEST_SOURCES)
  get_filename_component(test_name "${test_source}" NAME_WE)
  add_executable("${test_name}" "${test_source}")
  target_link_libraries("${test_name}" PRIVATE gss::gss GTest::gtest_main)
  target_compile_definitions("${test_name}" PRIVATE GSS_TEST_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  gss_apply_common_flags("${test_name}")
  add_test(NAME "${test_name}" COMMAND "${test_name}")
  if(VERIFICATION_SANITIZERS)
    set_tests_properties("${test_name}" PROPERTIES ENVIRONMENT
      "OMPI_MCA_memory=^patcher;ASAN_OPTIONS=detect_leaks=1:halt_on_error=1;UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1;LSAN_OPTIONS=suppressions=${VERIFICATION_SUPPRESSION_FILE}:print_suppressions=1;OMP_NUM_THREADS=1")
  endif()
endforeach()
