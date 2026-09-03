# Spikes: demoted code, off by default. See README.md for the rules.
#
# Deliberately minimal. These targets exist so the code is compiled often enough
# to notice when it rots, not so anything can depend on them; nothing here is
# installed and nothing here is registered with ctest.

message(STATUS "SCCD: building spikes (not installed, not tested)")

# The CUDA spikes are self-contained translation units. They are compiled into a
# separate object library rather than into sccd, so that a spike can never end up
# linked into the shipped library by accident.
if(SCCD_ENABLE_CUDA)
  file(GLOB SCCD_SPIKE_CUDA_SOURCES CONFIGURE_DEPENDS
       "${CMAKE_CURRENT_SOURCE_DIR}/src/cuda/*.cu")
  if(SCCD_SPIKE_CUDA_SOURCES)
    add_library(sccd_spikes_cuda OBJECT ${SCCD_SPIKE_CUDA_SOURCES})
    target_include_directories(sccd_spikes_cuda PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
      "${CMAKE_CURRENT_SOURCE_DIR}/src/cuda")
    target_link_libraries(sccd_spikes_cuda PRIVATE sccd)
    set_target_properties(sccd_spikes_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()

  # The demo that is the lower-bound kernel's only caller. It needs smesh like
  # the other demos do.
  if(SCCD_ENABLE_SMESH AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/demo/cuda/mesh_lower_bound_cuda.exe.cpp")
    add_executable(mesh_lower_bound_cuda "${CMAKE_CURRENT_SOURCE_DIR}/demo/cuda/mesh_lower_bound_cuda.exe.cpp")
    target_include_directories(mesh_lower_bound_cuda PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
      "${CMAKE_CURRENT_SOURCE_DIR}/src/cuda")
    target_link_libraries(mesh_lower_bound_cuda PRIVATE sccd)
    if(TARGET sccd_spikes_cuda)
      target_link_libraries(mesh_lower_bound_cuda PRIVATE sccd_spikes_cuda)
    endif()
  endif()
endif()

# Header-only spikes are compiled by a translation unit that does nothing but
# include them. Without this they would never be compiled at all, and a header
# nobody compiles is a header that quietly stops building.
file(GLOB SCCD_SPIKE_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp")
if(SCCD_SPIKE_HEADERS)
  set(_spike_tu "${CMAKE_CURRENT_BINARY_DIR}/spike_headers.cpp")
  set(_spike_tu_content "// Generated. Compiles every header-only spike so it cannot rot.\n")
  foreach(_spike_header IN LISTS SCCD_SPIKE_HEADERS)
    get_filename_component(_spike_header_name "${_spike_header}" NAME)
    string(APPEND _spike_tu_content "#include \"${_spike_header_name}\"\n")
  endforeach()
  file(GENERATE OUTPUT "${_spike_tu}" CONTENT "${_spike_tu_content}")
  add_library(sccd_spikes_headers OBJECT "${_spike_tu}")
  target_include_directories(sccd_spikes_headers PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
  target_link_libraries(sccd_spikes_headers PRIVATE sccd)
endif()

# hacks/ is not added here. It cannot configure -- it references a tests/
# directory that does not exist and needs scalable_ccd, Catch2 and libigl -- and
# wiring it up would mean either fetching those or letting the spike build fail,
# both of which defeat the point of the switch.
