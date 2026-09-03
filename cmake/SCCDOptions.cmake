option(SCCD_ENABLE_TBB "Use TBB for parallelization" OFF)
option(SCCD_ENABLE_OPENMP "Use OpenMP for parallelization" ON)
option(SCCD_ENABLE_CUDA "Enable CUDA" OFF)
# Left unset, CMake picks a default architecture for the detected toolkit,
# which on a recent toolkit is well below the target hardware. Register file
# size, shared-memory limits and DFMA scheduling all differ by architecture,
# so a kernel tuned against the default is tuned against the wrong machine.
# 90 is GH200 (Alps), the benchmark target.
set(SCCD_CUDA_ARCHITECTURES "90" CACHE STRING "CUDA architectures to build for")
option(SCCD_BUILD_SHARED_LIBS "Build SCCD as a shared library" ON)
option(SCCD_ENABLE_NATIVE_ARCH "Build for the host CPU (-march=native); required for the AVX2/AVX-512 kernels" ON)
option(SCCD_ENABLE_SPIKES "Build demoted code under spikes/ (not installed, not tested)" OFF)
option(SCCD_ENABLE_SMESH "Enable smesh and demos" OFF)
option(SCCD_ENABLE_TIGHT_INCLUSION "Enable Tight Inclusion" OFF)
option(SCCD_USE_VNARROW_PHASE_DEFAULT "Use the vectorized VF narrow phase by default" OFF)
option(SCCD_VNARROWPHASE_TI_COMPAT_DEFAULT "Correct vectorized VF outputs with TightInclusion by default" OFF)

if(SCCD_VNARROWPHASE_TI_COMPAT_DEFAULT AND NOT SCCD_ENABLE_TIGHT_INCLUSION)
  message(FATAL_ERROR "SCCD_VNARROWPHASE_TI_COMPAT_DEFAULT requires SCCD_ENABLE_TIGHT_INCLUSION=ON")
endif()

# Without an explicit build type CMake passes no optimization flags at all, which
# for this library means an -O0 build that looks like a normal one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
  message(STATUS "SCCD: CMAKE_BUILD_TYPE was unset, defaulting to Release")
endif()

if(SCCD_ENABLE_TBB AND SCCD_ENABLE_OPENMP)
  message(STATUS "SCCD: both TBB and OpenMP are enabled; TBB takes precedence in sparallel.hpp")
endif()
