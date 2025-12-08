cmake_minimum_required(VERSION 3.16)
project(sccd LANGUAGES CXX)

set(SCCD_VERSION_MAJOR 0)
set(SCCD_VERSION_MINOR 1)
set(SCCD_VERSION_PATCH 0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(SCCD_ENABLE_TBB "Use TBB for parallelization" ON)
option(SCCD_ENABLE_OPENMP "Use OpenMP for parallelization" OFF)
option(SCCD_ENABLE_HACKS "Do not use without Scalable-CCD library" OFF)

# SCCD library
file(GLOB_RECURSE SCCD_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")
file(GLOB_RECURSE SCCD_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(sccd SHARED ${SCCD_SOURCES})
target_sources(sccd PUBLIC ${SCCD_SOURCES};${SCCD_HEADERS})

if(SCCD_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(sccd PUBLIC OpenMP::OpenMP_CXX)
endif()

if(SCCD_ENABLE_TBB)
  find_package(TBB REQUIRED)
  if(TBB_FOUND)
    target_link_libraries(sccd PUBLIC TBB::tbb)
  endif()
endif()

if(SCCD_ENABLE_HACKS)
  file(GLOB_RECURSE SCCD_HACKS_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/hacks/*.hpp")
  file(GLOB_RECURSE SCCD_HACKS_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/hacks/*.cpp")
  target_sources(sccd PUBLIC ${SCCD_HACKS_SOURCES};${SCCD_HACKS_HEADERS})
  target_include_directories(sccd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/hacks>)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/sccd_config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/sccd_config.hpp")
target_include_directories(sccd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
target_include_directories(sccd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

# Installation
include(GNUInstallDirs)
install(TARGETS sccd
  EXPORT SCCD_TARGETS
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${SCCD_HEADERS} DESTINATION include)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python/sccd_py.py" DESTINATION python)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/sccd_config.hpp" DESTINATION include)

# install(FILES sccd.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sccd RENAME sccd.hpp)
