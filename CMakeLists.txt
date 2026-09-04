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

# NOTE: every SCCD_* option is declared in cmake/SCCDOptions.cmake, which is
# included above. Re-declaring them here would be a silent no-op.
# SCCD library
file(GLOB_RECURSE SCCD_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")
file(GLOB_RECURSE SCCD_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER SCCD_SOURCES EXCLUDE REGEX "/tests/")
if(NOT SCCD_ENABLE_SMESH)
  list(FILTER SCCD_SOURCES EXCLUDE REGEX "/integrations/smesh/")
  list(FILTER SCCD_HEADERS EXCLUDE REGEX "/integrations/smesh/")
endif()

if(SCCD_ENABLE_CUDA)
  file(GLOB_RECURSE SCCD_CUDA_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")
  file(GLOB_RECURSE SCCD_CUDA_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cuh")
  # A backup kept next to the translation unit it was copied from still ends in
  # .cu, so the glob compiles both and the linker sees every symbol twice.
  list(FILTER SCCD_CUDA_SOURCES EXCLUDE REGEX "\\.bkp\\.cu$")
  list(APPEND SCCD_SOURCES ${SCCD_CUDA_SOURCES})
  list(APPEND SCCD_HEADERS ${SCCD_CUDA_HEADERS})

  set(_sccd_library_include_dirs "")
  foreach(_file IN LISTS SCCD_CUDA_HEADERS)
    if(_file MATCHES "\\.(h|hh|hpp|hxx|inl|cuh)$")
      get_filename_component(_file_dir "${_file}" DIRECTORY)
      list(APPEND _sccd_library_include_dirs "${_file_dir}")
    endif()
  endforeach()

  list(REMOVE_DUPLICATES _sccd_library_include_dirs)
  set(_sccd_library_build_include_dirs "")
  foreach(_inc_dir IN LISTS _sccd_library_include_dirs)
      list(APPEND _sccd_library_build_include_dirs "$<BUILD_INTERFACE:${_inc_dir}>")
endforeach()

endif()

if(SCCD_BUILD_SHARED_LIBS)
  add_library(sccd SHARED ${SCCD_SOURCES})
else()
  add_library(sccd STATIC ${SCCD_SOURCES})
endif()
target_compile_features(sccd PUBLIC cxx_std_17)

# SCCD is almost entirely header-only, and the AVX2/AVX-512/NEON kernels in
# vaabb.hpp are guarded on __AVX2__/__AVX512F__/__ARM_NEON. Those macros only
# exist if the ISA flags are on the command line, so without this the hottest
# broadphase kernel silently degrades to the scalar fallback -- in consumer
# translation units as much as in ours. Hence PUBLIC, not PRIVATE.
if(SCCD_ENABLE_NATIVE_ARCH)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" SCCD_HAS_MARCH_NATIVE)
  if(SCCD_HAS_MARCH_NATIVE)
    target_compile_options(sccd PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
    message(STATUS "SCCD: building for the host CPU (-march=native)")
  else()
    check_cxx_compiler_flag("/arch:AVX2" SCCD_HAS_ARCH_AVX2)
    if(SCCD_HAS_ARCH_AVX2)
      target_compile_options(sccd PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/arch:AVX2>)
      message(STATUS "SCCD: building with /arch:AVX2")
    else()
      message(WARNING "SCCD: no host-arch flag available; the vectorized AABB kernels will fall back to scalar code")
    endif()
  endif()
else()
  message(WARNING
    "SCCD_ENABLE_NATIVE_ARCH=OFF: __AVX2__/__AVX512F__ will be undefined and the "
    "vectorized AABB overlap kernels will fall back to scalar code. Pass an explicit "
    "-march=/-mavx2 in CMAKE_CXX_FLAGS if you need a portable binary with SIMD.")
endif()

target_include_directories(sccd
    PUBLIC
        ${_sccd_library_build_include_dirs}
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

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

if(SCCD_ENABLE_CUDA)
  target_compile_definitions(sccd PUBLIC SCCD_ENABLE_CUDA)

  # The narrow-phase DFS kernels are register-bound: they carry a box, a set of
  # corner values and the per-query geometry live across the whole loop, and the
  # double instantiations sit near the point where nvcc starts spilling to local
  # memory. Spills there are invisible in a timing run but cost several times the
  # kernel's runtime, so make the register and spill counts printable.
  option(SCCD_CUDA_VERBOSE_PTXAS "Print per-kernel register and spill counts" OFF)
  if(SCCD_CUDA_VERBOSE_PTXAS)
    target_compile_options(sccd PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas=-v,-warn-spills>)
  endif()
endif()

# Demoted code. Off by default, never installed, never in the correctness gate;
# spikes/README.md states the rules and what is in there and why. hacks/ moved
# under it -- SCCD_ENABLE_HACKS is gone, since hacks/ cannot configure on its own
# and the spike switch is now the one knob for all of this.
if(SCCD_ENABLE_SPIKES)
  add_subdirectory(spikes)
endif()

if(SCCD_ENABLE_SMESH)
  target_include_directories(sccd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/integrations/smesh>)
  target_link_libraries(sccd PUBLIC smesh::smesh)
endif()

target_link_libraries(sccd PUBLIC ${SCCD_DEP_LIBRARIES})

include(GNUInstallDirs)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/core/sccd_config.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/sccd_config.hpp")
# Includes are written flat -- `#include "sccd_narrowphase.hpp"`, not
# `#include "narrowphase/sccd_narrowphase.hpp"` -- so every directory holding a
# header goes on the path, in the build tree and after install alike. The
# `sccd_` prefix plus one descriptive name per file keeps basenames unique,
# which is what makes flat includes safe.
set(_sccd_header_dirs core broadphase narrowphase cuda integrations/smesh)
set(_sccd_build_includes
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)
set(_sccd_install_includes $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
foreach(_d IN LISTS _sccd_header_dirs)
  list(APPEND _sccd_build_includes $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/${_d}>)
  list(APPEND _sccd_install_includes $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${_d}>)
endforeach()

target_include_directories(sccd
  PUBLIC
    ${_sccd_build_includes}
    ${_sccd_install_includes}
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

# The installed header set is listed, not globbed.
#
# A flat `install(FILES ${SCCD_HEADERS})` made every header under src/ public API
# by accident -- that is how broadphase_lb.hpp, which had no caller outside
# hacks/, came to be installed next to narrowphase.hpp. Listing them means adding
# a header is a decision someone makes rather than something that happens.
#
# Everything here is header-only templates, so the support headers are installed
# too: a consumer of narrowphase.hpp needs everything it includes.
set(SCCD_PUBLIC_HEADERS
  # Entry points.
  src/narrowphase/sccd_narrowphase.hpp
  src/narrowphase/sccd_narrowphase_quad.hpp
  src/broadphase/sccd_broadphase_sweep.hpp
  src/broadphase/sccd_broadphase_cell2d.hpp
  src/broadphase/sccd_broadphase_strategy.hpp
  # Configuration and shared vocabulary.
  src/core/sccd_base.hpp
  src/narrowphase/sccd_narrowphase_mode.hpp
  # Implementation headers, public because the entry points include them.
  src/narrowphase/sccd_objective.hpp
  src/narrowphase/sccd_rootfinder.hpp
  src/narrowphase/sccd_rootfinder_quad.hpp
  src/narrowphase/sccd_narrowphase_tight.hpp
  src/core/sccd_aabb.hpp
  src/core/sccd_math.hpp
  src/narrowphase/sccd_tolerance.hpp
  src/narrowphase/sccd_numerical_error.hpp
  src/core/sccd_parallel.hpp
)
if(SCCD_ENABLE_SMESH)
  list(APPEND SCCD_PUBLIC_HEADERS src/integrations/smesh/sccd_smesh_ccd.hpp)
endif()

# Fail configure if a header exists that nobody has classified. Without this the
# list silently drifts from the tree and we are back to guessing.
file(GLOB_RECURSE _sccd_all_headers RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")
list(FILTER _sccd_all_headers EXCLUDE REGEX "^src/tests/")
if(NOT SCCD_ENABLE_SMESH)
  list(FILTER _sccd_all_headers EXCLUDE REGEX "^src/integrations/smesh/")
endif()
foreach(_hdr IN LISTS _sccd_all_headers)
  if(NOT _hdr IN_LIST SCCD_PUBLIC_HEADERS)
    message(FATAL_ERROR
      "SCCD: ${_hdr} is not in SCCD_PUBLIC_HEADERS. Add it there if it is part of "
      "the public surface, or move it to spikes/ if it is not. See spikes/README.md.")
  endif()
endforeach()

foreach(_hdr IN LISTS SCCD_PUBLIC_HEADERS)
  get_filename_component(_hdr_dir "${_hdr}" DIRECTORY)
  string(REGEX REPLACE "^src/?" "" _hdr_dir "${_hdr_dir}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/${_hdr}"
          DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${_hdr_dir}")
endforeach()
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


# Demos that need nothing but the library itself build in every configuration --
# they are what a reader checks the API against, so gating them on an optional
# dependency would defeat the point.
add_executable(sccd_minimal "${CMAKE_CURRENT_SOURCE_DIR}/demo/sccd_minimal.exe.cpp")
target_link_libraries(sccd_minimal PRIVATE sccd)
add_test(NAME sccd_minimal COMMAND $<TARGET_FILE:sccd_minimal>)

if(SCCD_ENABLE_SMESH)
  file(GLOB DEMO_DRIVERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/demo/*.exe.cpp")
  list(FILTER DEMO_DRIVERS EXCLUDE REGEX "sccd_minimal\\.exe\\.cpp$")
  foreach(DEMO_DRIVER ${DEMO_DRIVERS})
    get_filename_component(DEMO_TARGET ${DEMO_DRIVER} NAME_WE)
    string(REGEX REPLACE "\\.exe$" "" DEMO_TARGET "${DEMO_TARGET}")
    add_executable(${DEMO_TARGET} ${DEMO_DRIVER})
    target_link_libraries(${DEMO_TARGET} PRIVATE sccd)
  endforeach()

  if(SCCD_ENABLE_CUDA)
    file(GLOB_RECURSE SCCD_CUDA_DEMOS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/demo/cuda/*.exe.cpp")
    foreach(SCCD_CUDA_DEMO ${SCCD_CUDA_DEMOS})
      get_filename_component(SCCD_CUDA_DEMO_TARGET ${SCCD_CUDA_DEMO} NAME_WE)
      string(REGEX REPLACE "\\.exe$" "" SCCD_CUDA_DEMO_TARGET "${SCCD_CUDA_DEMO_TARGET}")
      add_executable(${SCCD_CUDA_DEMO_TARGET} ${SCCD_CUDA_DEMO})
      target_link_libraries(${SCCD_CUDA_DEMO_TARGET} PRIVATE sccd)
    endforeach()
  endif()

  add_executable(sccd_bench "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/bench.exe.cpp")
  target_link_libraries(sccd_bench PRIVATE sccd)

  # Scaling study over surface refinement. Needs smesh for the mesh, refine and
  # skin operations, so it lives with the other smesh-gated targets.
  add_executable(sccd_refine_scaling "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/refine_scaling.exe.cpp")
  target_link_libraries(sccd_refine_scaling PRIVATE sccd)
endif()

# Accuracy oracle: compares every narrow-phase mode against TightInclusion on the
# CCD benchmark query sets. Needs TightInclusion but deliberately not smesh, so
# it can be built and run without the mesh stack.
if(SCCD_ENABLE_TIGHT_INCLUSION)
  add_executable(ti_oracle "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/ti_oracle.exe.cpp")
  target_link_libraries(ti_oracle PRIVATE sccd)
  if(SCCD_ENABLE_CUDA)
    # The oracle stays a plain C++ translation unit -- it only calls the host-side
    # narrow_phase_* entry points -- but it needs cuda_runtime.h for the buffer
    # management and libcudart to link.
    target_link_libraries(ti_oracle PRIVATE CUDA::cudart)
  endif()
endif()

include(SCCDTests)
