###
### Compile fortran module sources
###
set(SOURCES
    global_variables.f90
    opt_variables.f90
    performance_analyzer.f90
    nvtx.f90
    )
#if (USE_NVTX)
#   set(SOURCES ${SOURCES} nvtx.f90)
#endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
