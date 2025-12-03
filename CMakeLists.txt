cmake_minimum_required(VERSION 3.16)
project(sccd LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library
add_library(sccd SHARED sccd.cpp)
target_include_directories(sccd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Optional dependencies used by narrowphase.hpp; if found, link them
find_package(TBB REQUIRED)
if(TBB_FOUND)
  target_link_libraries(sccd PUBLIC TBB::tbb)
endif()

# find_package(Eigen3 REQUIRED)
# if(Eigen3_FOUND)
#   target_link_libraries(sccd PUBLIC Eigen3::Eigen)
# endif()

# Installation
include(GNUInstallDirs)
install(TARGETS sccd
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES sccd.cpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sccd RENAME sccd.hpp)

message(STATUS "Build the shared library with:")
message(STATUS "  mkdir -p build && cd build && cmake .. && cmake --build . --config Release")


