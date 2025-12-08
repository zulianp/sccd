Add these lines in the  scalable_ccd CMakeLists.txt

```cmake
set(SCCD_ENABLE_HACKS ON)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/sccd)
target_link_libraries(scalable_ccd PUBLIC sccd)
```

If TI is avilable in scalable_ccd also add this line

```cmake
target_compile_definitions(sccd PUBLIC "-DSCCD_ENABLE_TIGHT_INCLUSION=1")
```