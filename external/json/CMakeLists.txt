cmake_minimum_required(VERSION 3.16)
project(sccd_json_tools LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(
  simdjson
  GIT_REPOSITORY https://github.com/simdjson/simdjson.git
  GIT_TAG v3.13.0
)
FetchContent_MakeAvailable(simdjson)

add_executable(boxes_json_to_raw boxes_json_to_raw.cpp)
target_compile_features(boxes_json_to_raw PRIVATE cxx_std_17)
target_link_libraries(boxes_json_to_raw PRIVATE simdjson::simdjson)
