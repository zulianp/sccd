# ==============================================================================
# Tests
# ==============================================================================

add_executable(sccd_tests
  ${CMAKE_SOURCE_DIR}/tests/main.cpp
  test_sccd_narrow.cpp test_sccd_broad.cpp
  # Utils
  ${CMAKE_SOURCE_DIR}/tests/ground_truth.cpp
  ${CMAKE_SOURCE_DIR}/tests/ground_truth.hpp
  ${CMAKE_SOURCE_DIR}/tests/io.cpp
  ${CMAKE_SOURCE_DIR}/tests/io.hpp
)


# ==============================================================
# Required Libraries
# ==============================================================================
target_include_directories(sccd_tests PRIVATE ${CMAKE_SOURCE_DIR}/tests)
target_link_libraries(sccd_tests PRIVATE scalable_ccd::scalable_ccd)

# Catch2
include(catch2)
target_link_libraries(sccd_tests PRIVATE Catch2::Catch2)

# libigl
include(libigl)
target_link_libraries(sccd_tests PRIVATE igl::core)

# JSON
include(json)
target_link_libraries(sccd_tests PRIVATE nlohmann_json::nlohmann_json)

# TBB
include(onetbb)
target_link_libraries(sccd_tests PRIVATE TBB::tbb)

# Sample data
# include(sample_data)
# target_link_libraries(sccd_tests PRIVATE scalable_ccd::data)
set(SCALABLE_CCD_DATA_DIR "${CMAKE_SOURCE_DIR}/tests/data/") 
target_compile_definitions(sccd_tests PUBLIC  SCALABLE_CCD_DATA_DIR=\"${SCALABLE_CCD_DATA_DIR}\")

# Extra warnings (link last for highest priority)
include(scalable_ccd_warnings)
target_link_libraries(sccd_tests PRIVATE scalable_ccd::warnings)

# ==============================================================================
# Compiler options
# ==============================================================================

target_compile_definitions(sccd_tests PRIVATE CATCH_CONFIG_ENABLE_BENCHMARKING)

# ==============================================================================
# Register tests
# ==============================================================================

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(sccd_tests)
