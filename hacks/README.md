Modify the  scalable_ccd CMakeLists.txt

```cmake

set(SCCD_ENABLE_HACKS ON)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/sccd)

target_compile_definitions(sccd PUBLIC "-DSCCD_ENABLE_TIGHT_INCLUSION=1")

# Add an empty library and fill in the list of sources in `src/scalable_ccd/CMakeLists.txt`.
add_library(scalable_ccd)
add_library(scalable_ccd::scalable_ccd ALIAS scalable_ccd)
target_link_libraries(sccd PUBLIC Eigen3::Eigen)
target_include_directories(sccd PUBLIC ${SCALABLE_CCD_INCLUDE_DIR})

add_subdirectory("${SCALABLE_CCD_INCLUDE_DIR}/tight_inclusion")

target_link_libraries(scalable_ccd PUBLIC sccd)
```

If TI is avilable in scalable_ccd also add this line

```cmake
target_compile_definitions(sccd PUBLIC "-DSCCD_ENABLE_TIGHT_INCLUSION=1")
```