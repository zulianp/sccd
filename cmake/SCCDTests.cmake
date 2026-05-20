if(NOT SCCD_ENABLE_CUDA)
  return()
endif()

function(sccd_find_db_to_raw SCCD_DB_TO_RAW_OUT)
  if(SCCD_DB_TO_RAW_EXECUTABLE)
    set(${SCCD_DB_TO_RAW_OUT} "${SCCD_DB_TO_RAW_EXECUTABLE}" PARENT_SCOPE)
    return()
  endif()

  set(SCCD_SMESH_BIN_HINTS)
  if(TARGET smesh::smesh)
    get_target_property(SCCD_SMESH_IMPORTED_CONFIGS smesh::smesh IMPORTED_CONFIGURATIONS)
    if(SCCD_SMESH_IMPORTED_CONFIGS)
      foreach(SCCD_SMESH_IMPORTED_CONFIG IN LISTS SCCD_SMESH_IMPORTED_CONFIGS)
        get_target_property(SCCD_SMESH_LIBRARY smesh::smesh IMPORTED_LOCATION_${SCCD_SMESH_IMPORTED_CONFIG})
        if(SCCD_SMESH_LIBRARY)
          get_filename_component(SCCD_SMESH_LIBRARY_DIR "${SCCD_SMESH_LIBRARY}" DIRECTORY)
          get_filename_component(SCCD_SMESH_INSTALL_PREFIX "${SCCD_SMESH_LIBRARY_DIR}" DIRECTORY)
          list(APPEND SCCD_SMESH_BIN_HINTS "${SCCD_SMESH_INSTALL_PREFIX}/bin")
        endif()
      endforeach()
    endif()

    get_target_property(SCCD_SMESH_LIBRARY smesh::smesh IMPORTED_LOCATION)
    if(SCCD_SMESH_LIBRARY)
      get_filename_component(SCCD_SMESH_LIBRARY_DIR "${SCCD_SMESH_LIBRARY}" DIRECTORY)
      get_filename_component(SCCD_SMESH_INSTALL_PREFIX "${SCCD_SMESH_LIBRARY_DIR}" DIRECTORY)
      list(APPEND SCCD_SMESH_BIN_HINTS "${SCCD_SMESH_INSTALL_PREFIX}/bin")
    endif()
  endif()

  list(REMOVE_DUPLICATES SCCD_SMESH_BIN_HINTS)
  find_program(SCCD_DB_TO_RAW_EXECUTABLE db_to_raw HINTS ${SCCD_SMESH_BIN_HINTS} REQUIRED)
  set(${SCCD_DB_TO_RAW_OUT} "${SCCD_DB_TO_RAW_EXECUTABLE}" PARENT_SCOPE)
endfunction()

function(sccd_add_raw_mesh_test SCCD_TEST_NAME SCCD_TEST_TARGET SCCD_RAW_TARGET SCCD_RAW_DIR)
  set(SCCD_MESHES ${ARGN})
  if(NOT SCCD_MESHES)
    return()
  endif()

  foreach(SCCD_MESH IN LISTS SCCD_MESHES)
    if(NOT EXISTS "${SCCD_MESH}")
      return()
    endif()
  endforeach()

  sccd_find_db_to_raw(SCCD_DB_TO_RAW)
  set(SCCD_RAW_STAMP "${SCCD_RAW_DIR}/converted.stamp")
  set(SCCD_RAW_OUTPUTS)
  set(SCCD_RAW_COMMANDS)
  foreach(SCCD_MESH IN LISTS SCCD_MESHES)
    get_filename_component(SCCD_MESH_NAME "${SCCD_MESH}" NAME_WE)
    set(SCCD_RAW_OUTPUT "${SCCD_RAW_DIR}/${SCCD_MESH_NAME}")
    list(APPEND SCCD_RAW_OUTPUTS "${SCCD_RAW_OUTPUT}")
    list(APPEND SCCD_RAW_COMMANDS COMMAND "${SCCD_DB_TO_RAW}" "${SCCD_MESH}" "${SCCD_RAW_OUTPUT}")
  endforeach()

  add_custom_command(
    OUTPUT "${SCCD_RAW_STAMP}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${SCCD_RAW_DIR}"
    ${SCCD_RAW_COMMANDS}
    COMMAND ${CMAKE_COMMAND} -E touch "${SCCD_RAW_STAMP}"
    DEPENDS ${SCCD_MESHES}
    VERBATIM)
  add_custom_target(${SCCD_RAW_TARGET} DEPENDS "${SCCD_RAW_STAMP}")
  add_dependencies(${SCCD_TEST_TARGET} ${SCCD_RAW_TARGET})

  add_test(
    NAME ${SCCD_TEST_NAME}
    COMMAND $<TARGET_FILE:${SCCD_TEST_TARGET}> ${SCCD_RAW_OUTPUTS})
endfunction()

enable_testing()

file(GLOB_RECURSE SCCD_CUDA_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/tests/cuda/*.exe.cpp")
foreach(SCCD_CUDA_TEST IN LISTS SCCD_CUDA_TESTS)
  get_filename_component(SCCD_CUDA_TEST_TARGET "${SCCD_CUDA_TEST}" NAME_WE)
  string(REGEX REPLACE "\\.exe$" "" SCCD_CUDA_TEST_TARGET "${SCCD_CUDA_TEST_TARGET}")
  add_executable(${SCCD_CUDA_TEST_TARGET} "${SCCD_CUDA_TEST}")
  target_link_libraries(${SCCD_CUDA_TEST_TARGET} PRIVATE sccd)
  add_test(NAME ${SCCD_CUDA_TEST_TARGET} COMMAND $<TARGET_FILE:${SCCD_CUDA_TEST_TARGET}>)
endforeach()

if(TARGET mesh_sccd_cuda_test)
  sccd_add_raw_mesh_test(
    sccd_cuda_cpu_gpu_parity
    mesh_sccd_cuda_test
    sccd_cuda_cpu_gpu_parity_meshes
    "${CMAKE_CURRENT_BINARY_DIR}/tmp/sccd_cuda_cpu_gpu_parity"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/n-body-simulation/frames/balls16_18.ply"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/n-body-simulation/frames/balls16_19.ply")
endif()
