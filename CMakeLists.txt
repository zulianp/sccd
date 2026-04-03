cmake_minimum_required(VERSION 3.16)
project(sccd LANGUAGES C CXX)

set(SCCD_VERSION_MAJOR 0)
set(SCCD_VERSION_MINOR 1)
set(SCCD_VERSION_PATCH 0)
set(SCCD_VERSION "${SCCD_VERSION_MAJOR}.${SCCD_VERSION_MINOR}.${SCCD_VERSION_PATCH}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(SCCDOptions)
include(SCCDDependencies)

option(SCCD_ENABLE_TBB "Use TBB for parallelization" ON)
option(SCCD_ENABLE_OPENMP "Use OpenMP for parallelization" OFF)
option(SCCD_ENABLE_HACKS "Do not use without Scalable-CCD library" OFF)
option(SCCD_ENABLE_SMESH "Enable smesh and demos" OFF)

# SCCD library
file(GLOB_RECURSE SCCD_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")
file(GLOB_RECURSE SCCD_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(sccd SHARED ${SCCD_SOURCES})
add_library(SCCD::sccd ALIAS sccd)
set_target_properties(sccd PROPERTIES
  EXPORT_NAME sccd
  VERSION "${SCCD_VERSION}"
  SOVERSION "${SCCD_VERSION_MAJOR}"
)

if(TARGET OpenMP::OpenMP_CXX)
  target_link_libraries(sccd PUBLIC OpenMP::OpenMP_CXX)
endif()

if(TARGET TBB::tbb)
  target_link_libraries(sccd PUBLIC TBB::tbb)
endif()

if(SCCD_ENABLE_HACKS)
  add_subdirectory(hacks)
  # file(GLOB_RECURSE SCCD_HACKS_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/hacks/*.hpp")
  # file(GLOB_RECURSE SCCD_HACKS_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/hacks/*.cpp")
  # target_sources(sccd PUBLIC ${SCCD_HACKS_SOURCES};${SCCD_HACKS_HEADERS})
  # target_include_directories(sccd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/hacks>)
endif()

if(SCCD_ENABLE_SMESH)
  target_link_libraries(sccd PUBLIC smesh::smesh)
endif()

include(GNUInstallDirs)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/sccd_config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/sccd_config.hpp")
target_include_directories(sccd
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Installation
include(CMakePackageConfigHelpers)

set(SCCD_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/SCCD")

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SCCDConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/SCCDConfig.cmake"
  INSTALL_DESTINATION "${SCCD_INSTALL_CMAKEDIR}"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/SCCDConfigVersion.cmake"
  VERSION "${SCCD_VERSION}"
  COMPATIBILITY SameMajorVersion
)

install(TARGETS sccd
  EXPORT SCCD_TARGETS
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES ${SCCD_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python/sccd_py.py" DESTINATION python)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/sccd_config.hpp" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT SCCD_TARGETS
  FILE SCCDTargets.cmake
  NAMESPACE SCCD::
  DESTINATION "${SCCD_INSTALL_CMAKEDIR}"
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SCCDConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SCCDConfigVersion.cmake"
  DESTINATION "${SCCD_INSTALL_CMAKEDIR}"
)

# install(FILES sccd.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/sccd RENAME sccd.hpp)


if(SCCD_ENABLE_SMESH)
  file(GLOB DEMO_DRIVERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/demo/*.exe.cpp")
  foreach(DEMO_DRIVER ${DEMO_DRIVERS})
    get_filename_component(DEMO_TARGET ${DEMO_DRIVER} NAME_WE)
    string(REGEX REPLACE "\\.exe$" "" DEMO_TARGET "${DEMO_TARGET}")
    add_executable(${DEMO_TARGET} ${DEMO_DRIVER})
    target_link_libraries(${DEMO_TARGET} PRIVATE sccd)
  endforeach()
endif()
