if(SCCD_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if(SCCD_ENABLE_TBB)
  find_package(TBB REQUIRED)
endif()

if(SCCD_ENABLE_CUDA)
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
  find_package(smesh REQUIRED)
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
