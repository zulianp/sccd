# Competitor harnesses: third-party CCD run over the same benchmark.
#
# Off by default. The shipped library builds with no external dependencies, and
# these pull Eigen, oneTBB and spdlog through the competitor's own CMake; keeping
# them behind a switch means the default configuration, and the correctness gate,
# stay free of all of it.
#
# benchmark/competitors/README.md states the rules a comparison here follows.

include(FetchContent)

# The harnesses read the prepared datasets through smesh, exactly as
# benchmark/bench.exe.cpp does, so that both sides see the same geometry from the
# same files.
if(NOT SCCD_ENABLE_SMESH)
  message(WARNING
    "SCCD_ENABLE_COMPETITORS=ON needs SCCD_ENABLE_SMESH=ON to read the benchmark "
    "meshes; skipping the competitor harnesses.")
  return()
endif()

# ---------------------------------------------------------------------------
# Scalable CCD
# ---------------------------------------------------------------------------

# Pinned to a commit, not to a branch: a comparison that cannot be reproduced a
# year from now is not worth publishing.
#
# Deliberately not to a tag. The newest release, v0.1.1, is an earlier library
# under a different name -- `project(CCDGPU LANGUAGES CXX CUDA)`, a libigl
# dependency, CUDA not optional, and none of the `scalable_ccd::` interface the
# paper describes. The broad and narrow phases only became separately callable
# afterwards, and calling them apart is what lets the phase split here match the
# article's. A commit hash pins just as firmly as a tag.
set(SCCD_SCALABLE_CCD_TAG "ccf4bc19aba432b9d5abdd41ae55d15ea166a0a5"
    CACHE STRING "Scalable CCD commit to compare against")

# Its CUDA support follows ours, so a CUDA build of SCCD is compared against a
# CUDA build of Scalable CCD and a host build against a host build.
set(SCALABLE_CCD_WITH_CUDA ${SCCD_ENABLE_CUDA} CACHE BOOL "" FORCE)
# Off is the library as its authors ship it: the narrow phase prunes every box
# against the step's running earliest time of impact and returns that one
# number. On, each query keeps its own bound and reports its own time of impact,
# which is what per-query accuracy needs and is strictly more work. So a build
# with it off is the cost measurement, a build with it on is the per-query
# accuracy measurement, and run_comparison.sh builds both.
option(SCCD_SCALABLE_CCD_TOI_PER_QUERY
  "Build Scalable CCD with per-query times of impact (accuracy runs; not its cost)" OFF)
set(SCALABLE_CCD_TOI_PER_QUERY ${SCCD_SCALABLE_CCD_TOI_PER_QUERY} CACHE BOOL "" FORCE)
set(SCALABLE_CCD_WITH_PROFILER OFF CACHE BOOL "" FORCE)
set(SCALABLE_CCD_BUILD_TESTS OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  scalable_ccd
  GIT_REPOSITORY https://github.com/continuous-collision-detection/scalable-ccd.git
  GIT_TAG ${SCCD_SCALABLE_CCD_TAG}
)
FetchContent_MakeAvailable(scalable_ccd)

# The same flags SCCD is built with, for the same reason cmake/SCCDDependencies
# gives TightInclusion them: a competitor left at the toolchain's defaults while
# SCCD is built for the host CPU measures the flags, not the algorithms.
if(TARGET scalable_ccd)
  separate_arguments(SCCD_COMPETITOR_CXX_FLAGS NATIVE_COMMAND "${CMAKE_CXX_FLAGS}")
  separate_arguments(SCCD_COMPETITOR_CXX_RELEASE_FLAGS NATIVE_COMMAND "${CMAKE_CXX_FLAGS_RELEASE}")
  target_compile_options(scalable_ccd PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SCCD_COMPETITOR_CXX_FLAGS}>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:${SCCD_COMPETITOR_CXX_RELEASE_FLAGS}>
  )
endif()

add_executable(scalable_ccd_bench
  "${CMAKE_CURRENT_SOURCE_DIR}/scalable_ccd/scalable_ccd_bench.exe.cpp")
target_include_directories(scalable_ccd_bench PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/common")
target_link_libraries(scalable_ccd_bench PRIVATE sccd scalable_ccd::scalable_ccd)

if(SCCD_ENABLE_CUDA)
  # The harness calls the device broad and narrow phases directly so it can time
  # them apart, which means it is a CUDA translation unit.
  set_source_files_properties(
    "${CMAKE_CURRENT_SOURCE_DIR}/scalable_ccd/scalable_ccd_bench.exe.cpp"
    PROPERTIES LANGUAGE CUDA)
  target_compile_definitions(scalable_ccd_bench PRIVATE SCCD_COMPETITOR_WITH_CUDA)
  set_target_properties(scalable_ccd_bench PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

# ---------------------------------------------------------------------------
# Additive CCD, from the IPC toolkit
# ---------------------------------------------------------------------------
#
# The method rather than a transcription of it: the rule here is that a
# competitor is the code its authors published, and the toolkit is where the
# preprint's implementation lives.
set(SCCD_IPC_TOOLKIT_TAG "v1.6.0" CACHE STRING "IPC toolkit version to compare against")

set(IPC_TOOLKIT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_CUDA OFF CACHE BOOL "" FORCE)
# Additive CCD needs none of these, and each one is a dependency to resolve
# against the ones Scalable CCD already fetched into the same tree.
set(IPC_TOOLKIT_WITH_SIMD OFF CACHE BOOL "" FORCE)
set(IPC_TOOLKIT_WITH_ABSEIL OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
  ipc_toolkit
  GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
  GIT_TAG ${SCCD_IPC_TOOLKIT_TAG}
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(ipc_toolkit)

if(TARGET ipc_toolkit)
  target_compile_options(ipc_toolkit PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${SCCD_COMPETITOR_CXX_FLAGS}>
    $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:${SCCD_COMPETITOR_CXX_RELEASE_FLAGS}>
  )
endif()

# A narrow phase over the curated queries, so it needs neither CUDA nor a mesh.
add_executable(accd_bench "${CMAKE_CURRENT_SOURCE_DIR}/accd/accd_bench.exe.cpp")
target_include_directories(accd_bench PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/common")
target_link_libraries(accd_bench PRIVATE sccd ipc::toolkit)

# Additive CCD answers one pair at a time, so the harness drives it over the
# candidates with the same `tbb::parallel_for` the toolkit's own stepsize query
# uses. oneTBB arrives through the toolkit, but its headers are not guaranteed to
# come with it -- they do on one of the two machines here and not on the other --
# so the dependency is named rather than inherited.
if(TARGET TBB::tbb)
  target_link_libraries(accd_bench PRIVATE TBB::tbb)
elseif(TARGET tbb)
  target_link_libraries(accd_bench PRIVATE tbb)
else()
  message(FATAL_ERROR "accd_bench needs oneTBB, which the IPC toolkit should have provided")
endif()

# The two binaries have to agree on the column list, or their rows cannot be
# concatenated -- which is the only reason to run the competitor at all. Each
# prints its own schema with --header; this compares them.
if(TARGET sccd_bench)
  foreach(_competitor scalable_ccd_bench accd_bench)
    add_test(NAME competitor_schema_${_competitor}
      COMMAND ${CMAKE_COMMAND}
        -DSCCD_BENCH=$<TARGET_FILE:sccd_bench>
        -DCOMPETITOR_BENCH=$<TARGET_FILE:${_competitor}>
        -P "${CMAKE_CURRENT_SOURCE_DIR}/compare_schema.cmake")
  endforeach()
endif()
