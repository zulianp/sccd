# The HIP platform must be known before enable_language(HIP): chipStar (HIP on SPIR-V,
# Intel GPUs) runs warps of 32 lanes where the AMD default is 64, and builds for no
# architecture list, and both are fixed before the language is enabled. Given, it is
# CMAKE_HIP_PLATFORM. Not given, enable_language takes it from `hipconfig --platform`
# (CMake's CMakeDetermineHIPCompiler.cmake, spirv from CMake 4.3), which ROCm's hipconfig
# answers from the HIP_PLATFORM environment variable when that is set. If either says
# spirv the configure stops here, before the AMD defaults are taken, and asks for the
# platform. SFEM does the same.
if(SCCD_ENABLE_HIP AND NOT CMAKE_HIP_PLATFORM)
    set(_sccd_hipconfig_platform "")
    execute_process(COMMAND hipconfig --platform
        OUTPUT_VARIABLE _sccd_hipconfig_output OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _sccd_hipconfig_result
        ERROR_QUIET)
    if(_sccd_hipconfig_result EQUAL 0)
        set(_sccd_hipconfig_platform "${_sccd_hipconfig_output}")
    endif()
    if("$ENV{HIP_PLATFORM}" STREQUAL "spirv" OR _sccd_hipconfig_platform STREQUAL "spirv")
        message(FATAL_ERROR "The HIP platform is chipStar's SPIR-V (HIP_PLATFORM: '$ENV{HIP_PLATFORM}', "
                            "hipconfig --platform: '${_sccd_hipconfig_platform}'), but CMAKE_HIP_PLATFORM was "
                            "not given: configure with -DCMAKE_HIP_PLATFORM=spirv")
    endif()
endif()
