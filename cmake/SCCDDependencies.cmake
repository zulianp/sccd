if(SCCD_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if(SCCD_ENABLE_TBB)
  find_package(TBB REQUIRED)
endif()

if(SCCD_ENABLE_HIP AND NOT SCCD_ENABLE_CUDA)
    message(FATAL_ERROR "SCCD_ENABLE_HIP=ON selects the backend of the device path, "
                        "which SCCD_ENABLE_CUDA=ON turns on: set both.")
endif()

if(SCCD_ENABLE_CUDA AND SCCD_ENABLE_HIP)
    # The device sources are compiled as HIP. Like the CUDA architectures below,
    # this must precede enable_language.
    if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
      set(CMAKE_HIP_ARCHITECTURES "${SCCD_HIP_ARCHITECTURES}")
    endif()

    enable_language(HIP)
    message(STATUS "SCCD: building the device path with HIP for ${CMAKE_HIP_ARCHITECTURES}, "
                   "wavefront size ${SCCD_HIP_WAVEFRONT_SIZE}")

    if(NOT DEFINED CMAKE_HIP_STANDARD)
      set(CMAKE_HIP_STANDARD 17)
      set(CMAKE_HIP_STANDARD_REQUIRED ON)
    endif()

    # Executables made of C++ sources are linked by the C++ compiler. A static
    # library with HIP objects -- smesh's, imported -- would otherwise hand them to
    # the HIP compiler, whose linker produces position-independent executables
    # from objects the host compiler built without -fPIE, and the link fails.
    set(CMAKE_HIP_LINKER_PREFERENCE_PROPAGATES OFF)

    find_package(hip REQUIRED)
    find_package(hipcub REQUIRED)
    find_package(rocthrust REQUIRED)

    # hip::host carries the runtime and the platform define that host code
    # including sccd_gpu_runtime.h needs, so it is public. hipCUB and rocThrust
    # are header-only and used inside the .cu files only; their targets are not
    # linked, because through rocprim they bring hip::device, which compiles every
    # C++ source of the library as HIP. CMakeLists.txt adds their include
    # directories to the HIP sources instead.
    list(APPEND SCCD_DEP_LIBRARIES hip::host)
elseif(SCCD_ENABLE_CUDA)
    # Must precede enable_language(CUDA): CMake fixes the default architecture
    # while enabling the language, and a later assignment is ignored for the
    # compiler test.
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES "${SCCD_CUDA_ARCHITECTURES}")
    endif()

    enable_language(CUDA)
    message(STATUS "SCCD: building CUDA for architectures ${CMAKE_CUDA_ARCHITECTURES}")

    if(NOT DEFINED CMAKE_CUDA_STANDARD)
      set(CMAKE_CUDA_STANDARD 17)
      set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()

    set(SCCD_DEP_INCLUDES "${SCCD_DEP_INCLUDES};${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
    #set(SCCD_DEP_LIBRARIES "${SCCD_DEP_LIBRARIES};${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}")

    include(CheckLanguage)
    check_language(CUDA)

    find_package(CUDAToolkit REQUIRED)

    list(APPEND SCCD_DEP_LIBRARIES "CUDA::cudart")

    set(_SCCD_CUDA_MODULES "CUDA::cusparse;CUDA::cublas;CUDA::nvToolsExt")

    set(SCCD_CUDA_MATH_LIBS_FOUND FALSE)

    foreach(CUDA_MODULE ${_SCCD_CUDA_MODULES})
        if(TARGET ${CUDA_MODULE})
            list(APPEND SCCD_DEP_LIBRARIES "${CUDA_MODULE}")
            set(SCCD_CUDA_MATH_LIBS_FOUND TRUE)
        else()
            message(WARNING "[Warning] CUDAToolkit does not have module ${CUDA_MODULE} in a standard location!")
            
        endif()
    endforeach()

    if(NOT SCCD_CUDA_MATH_LIBS_FOUND)
        message("Trying with: CRAY_CUDATOOLKIT_POST_LINK_OPTS=$ENV{CRAY_CUDATOOLKIT_POST_LINK_OPTS}")
        message("Trying with: CRAY_CUDATOOLKIT_INCLUDE_OPTS=$ENV{CRAY_CUDATOOLKIT_INCLUDE_OPTS}")
        list(APPEND SCCD_DEP_LIBRARIES "$ENV{CRAY_CUDATOOLKIT_POST_LINK_OPTS} -lcublas -lcusparse")
        include_directories($ENV{CRAY_CUDATOOLKIT_INCLUDE_OPTS})
    endif()

    # https://github.com/NVIDIA/thrust/blob/main/thrust/cmake/README.md
    #
    # PRIVATE, not PUBLIC: thrust is used inside src/cuda/sccd_narrowphase.cu and
    # by no installed header. Exporting it puts a target name in
    # SCCDTargets.cmake that a consumer has no way to define, and CMake then
    # falls back to a literal -lThrust, which does not exist.
    find_package(Thrust CONFIG)
    if(Thrust_FOUND)
        thrust_create_target(Thrust)
        list(APPEND SCCD_DEP_PRIVATE_LIBRARIES Thrust)
    else()
        message(WARNING "Thrust not found!")
    endif()
endif()

if(SCCD_ENABLE_SMESH)
  if(NOT TARGET smesh::smesh)
    find_package(smesh REQUIRED)
  else()
    # smesh built in this tree: its options are the backend.
    set(smesh_ENABLE_CUDA ${SMESH_ENABLE_CUDA})
    set(smesh_ENABLE_HIP ${SMESH_ENABLE_HIP})
  endif()

  # SCCD's device mesh path takes its points and connectivity from smesh's device
  # buffers, so smesh must have been built for the same device backend. An
  # install older than these variables says nothing and is not checked for CUDA;
  # it cannot be a HIP build.
  if(SCCD_ENABLE_CUDA AND DEFINED smesh_ENABLE_CUDA AND NOT smesh_ENABLE_CUDA)
    message(FATAL_ERROR "SCCD_ENABLE_CUDA=ON needs smesh built with SMESH_ENABLE_CUDA=ON: "
                        "without it the device mesh path receives host memory.")
  endif()
  if(SCCD_ENABLE_HIP AND NOT smesh_ENABLE_HIP)
    message(FATAL_ERROR "SCCD_ENABLE_HIP=ON needs smesh built with "
                        "SMESH_ENABLE_CUDA=ON and SMESH_ENABLE_HIP=ON.")
  endif()
  if(smesh_ENABLE_HIP AND NOT SCCD_ENABLE_HIP)
    message(FATAL_ERROR "smesh was built with HIP: configure SCCD with "
                        "SCCD_ENABLE_CUDA=ON and SCCD_ENABLE_HIP=ON.")
  endif()
endif()


if(SCCD_ENABLE_TIGHT_INCLUSION)
  include(FetchContent)

  set(TIGHT_INCLUSION_WITH_RATIONAL OFF CACHE BOOL "" FORCE)
  set(TIGHT_INCLUSION_WITH_TIMER OFF CACHE BOOL "" FORCE)
  set(TIGHT_INCLUSION_LIMIT_QUEUE_SIZE OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(
    tight_inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG v1.0.6
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(tight_inclusion)

  separate_arguments(SCCD_TIGHT_INCLUSION_CXX_FLAGS NATIVE_COMMAND "${CMAKE_CXX_FLAGS}")
  separate_arguments(SCCD_TIGHT_INCLUSION_CXX_RELEASE_FLAGS NATIVE_COMMAND "${CMAKE_CXX_FLAGS_RELEASE}")
  target_compile_options(tight_inclusion PRIVATE
    ${SCCD_TIGHT_INCLUSION_CXX_FLAGS}
    "$<$<CONFIG:Release>:${SCCD_TIGHT_INCLUSION_CXX_RELEASE_FLAGS}>"
  )

  list(APPEND SCCD_DEP_LIBRARIES "$<BUILD_INTERFACE:tight_inclusion::tight_inclusion>")

endif()
