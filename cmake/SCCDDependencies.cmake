if(SCCD_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if(SCCD_ENABLE_TBB)
  find_package(TBB REQUIRED)
endif()

if(SCCD_ENABLE_CUDA)
    enable_language(CUDA)

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

    #https://github.com/NVIDIA/thrust/blob/main/thrust/cmake/README.md
    find_package(Thrust CONFIG)
    if(Thrust_FOUND)
        thrust_create_target(Thrust)
        list(APPEND SCCD_DEP_LIBRARIES Thrust)
    else()
        message(WARNING "Thrust not found!")
    endif()
endif()

if(SCCD_ENABLE_SMESH)
  find_package(smesh REQUIRED)
endif()
