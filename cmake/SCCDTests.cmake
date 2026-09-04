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

  # Missing data used to mean this test quietly did not exist, which on a fresh
  # clone is always. Say so instead.
  foreach(SCCD_MESH IN LISTS SCCD_MESHES)
    if(NOT EXISTS "${SCCD_MESH}")
      message(STATUS "SCCD: not registering ${SCCD_TEST_NAME} -- missing ${SCCD_MESH}")
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

# Tests that need TightInclusion but not smesh. Listed by name so they can be
# excluded from the smesh glob below; otherwise enabling both options registers
# each of them twice and configuration fails on the duplicate target.
set(SCCD_TI_ONLY_TESTS
  sccd_numerical_error_ti_parity_test
  sccd_tolerance_ti_parity_test)

if(SCCD_ENABLE_TIGHT_INCLUSION)
  foreach(SCCD_TI_TEST IN LISTS SCCD_TI_ONLY_TESTS)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/${SCCD_TI_TEST}.exe.cpp")
      add_executable(${SCCD_TI_TEST} "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/${SCCD_TI_TEST}.exe.cpp")
      target_link_libraries(${SCCD_TI_TEST} PRIVATE sccd)
      add_test(NAME ${SCCD_TI_TEST} COMMAND $<TARGET_FILE:${SCCD_TI_TEST}>)
    endif()
  endforeach()
endif()

# Tests that need neither smesh nor TightInclusion, so they run in any build.
#
set(SCCD_STANDALONE_TESTS
  sccd_broadphase_cell2d_test
  sccd_rootfinder_quad_test
  sccd_c_abi_test
  sccd_aabb_test)
foreach(SCCD_TEST IN LISTS SCCD_STANDALONE_TESTS)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/${SCCD_TEST}.exe.cpp")
    add_executable(${SCCD_TEST} "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/${SCCD_TEST}.exe.cpp")
    target_link_libraries(${SCCD_TEST} PRIVATE sccd)
    add_test(NAME ${SCCD_TEST} COMMAND $<TARGET_FILE:${SCCD_TEST}>)
  endif()
endforeach()

if(SCCD_ENABLE_SMESH)
  file(GLOB SCCD_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/*.exe.cpp")
  foreach(SCCD_TEST IN LISTS SCCD_TESTS)
    get_filename_component(SCCD_TEST_TARGET "${SCCD_TEST}" NAME_WE)
    string(REGEX REPLACE "\\.exe$" "" SCCD_TEST_TARGET "${SCCD_TEST_TARGET}")
    # Already registered above, and it must stay gated on TightInclusion.
    if(SCCD_TEST_TARGET IN_LIST SCCD_TI_ONLY_TESTS OR SCCD_TEST_TARGET IN_LIST SCCD_STANDALONE_TESTS)
      continue()
    endif()
    add_executable(${SCCD_TEST_TARGET} "${SCCD_TEST}")
    target_link_libraries(${SCCD_TEST_TARGET} PRIVATE sccd)
    add_test(NAME ${SCCD_TEST_TARGET} COMMAND $<TARGET_FILE:${SCCD_TEST_TARGET}>)
  endforeach()
endif()

# The accuracy gate, as an actual test.
#
# ti_oracle is what decides whether a narrow-phase kernel is conservative -- it
# exits non-zero when a mode misses a collision or reports a time of impact later
# than the dataset's exact root -- and it was not a ctest, so the gate only ran
# when somebody remembered to run it.
#
# It needs the benchmark datasets, which are several gigabytes and not in the
# repository, so each scene registers only if its queries are actually present.
# --max-files bounds the runtime: the gate is meant to run with the suite, and a
# full trajectory per scene does not belong there. Run ti_oracle by hand for the
# exhaustive pass.
if(TARGET ti_oracle)
  set(SCCD_ORACLE_SCENES cloth-ball armadillo-rollers cloth-funnel)
  set(SCCD_ORACLE_MAX_FILES 8 CACHE STRING "Query files per scene in the ti_oracle ctest")
  foreach(SCCD_ORACLE_SCENE IN LISTS SCCD_ORACLE_SCENES)
    set(SCCD_ORACLE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/${SCCD_ORACLE_SCENE}")
    if(IS_DIRECTORY "${SCCD_ORACLE_DIR}/queries")
      add_test(
        NAME oracle_${SCCD_ORACLE_SCENE}
        COMMAND $<TARGET_FILE:ti_oracle> "${SCCD_ORACLE_DIR}"
                --max-files ${SCCD_ORACLE_MAX_FILES})
    else()
      message(STATUS
        "SCCD: not registering oracle_${SCCD_ORACLE_SCENE} -- no queries under ${SCCD_ORACLE_DIR}")
    endif()
  endforeach()
endif()

if(NOT SCCD_ENABLE_CUDA)
  return()
endif()

# CUDA tests that also need smesh. Registering these on SCCD_ENABLE_CUDA alone
# meant a CUDA build without smesh failed to compile them.
set(SCCD_CUDA_SMESH_TESTS sccd_mesh_cuda_test)

# CUDA tests that take arguments, and so must not be registered bare.
#
# sccd_mesh_cuda_test needs two mesh paths and aborts with its usage message
# without them. It is registered properly further down by sccd_add_raw_mesh_test,
# which supplies the meshes and skips when the data is absent -- so the bare
# registration below could only ever abort, and did: a CUDA build without the
# n-body dataset failed ctest on it every time. Neither local configuration
# showed it, because the machine with the data has no CUDA and the machine with
# CUDA has no data.
set(SCCD_CUDA_ARG_TESTS sccd_mesh_cuda_test)

file(GLOB_RECURSE SCCD_CUDA_TESTS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/cuda/*.exe.cpp")
foreach(SCCD_CUDA_TEST IN LISTS SCCD_CUDA_TESTS)
  get_filename_component(SCCD_CUDA_TEST_TARGET "${SCCD_CUDA_TEST}" NAME_WE)
  string(REGEX REPLACE "\\.exe$" "" SCCD_CUDA_TEST_TARGET "${SCCD_CUDA_TEST_TARGET}")
  if(SCCD_CUDA_TEST_TARGET IN_LIST SCCD_CUDA_SMESH_TESTS AND NOT SCCD_ENABLE_SMESH)
    message(STATUS "SCCD: skipping ${SCCD_CUDA_TEST_TARGET} (needs SCCD_ENABLE_SMESH)")
    continue()
  endif()
  add_executable(${SCCD_CUDA_TEST_TARGET} "${SCCD_CUDA_TEST}")
  target_link_libraries(${SCCD_CUDA_TEST_TARGET} PRIVATE sccd)
  # Built either way; only registered here if it can run without arguments.
  if(NOT SCCD_CUDA_TEST_TARGET IN_LIST SCCD_CUDA_ARG_TESTS)
    add_test(NAME ${SCCD_CUDA_TEST_TARGET} COMMAND $<TARGET_FILE:${SCCD_CUDA_TEST_TARGET}>)
  endif()
endforeach()

if(TARGET sccd_mesh_cuda_test)
  sccd_add_raw_mesh_test(
    sccd_cuda_cpu_gpu_parity
    sccd_mesh_cuda_test
    sccd_cuda_cpu_gpu_parity_meshes
    "${CMAKE_CURRENT_BINARY_DIR}/tmp/sccd_cuda_cpu_gpu_parity"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/n-body-simulation/frames/balls16_18.ply"
    "${CMAKE_CURRENT_SOURCE_DIR}/data/n-body-simulation/frames/balls16_19.ply")
endif()
