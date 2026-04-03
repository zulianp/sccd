if(SCCD_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if(SCCD_ENABLE_TBB)
  find_package(TBB REQUIRED)
endif()

if(SCCD_ENABLE_SMESH)
  find_package(smesh REQUIRED)
endif()
