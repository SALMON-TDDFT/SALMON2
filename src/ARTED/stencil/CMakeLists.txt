###
### Compile stencil computation sources
###

### Soruces for current and total_energy
set(SOURCES_EXPLICIT_VEC
    C/${SIMD_SET}/current.c
    C/${SIMD_SET}/total_energy.c
    )
set(SOURCES_COMPILER_VEC
    F90/total_energy.f90
    F90/current.f90
    )

###
if (USE_OPT_EXPLICIT_VECTORIZATION)
  set(SOURCES ${SOURCES_EXPLICIT_VEC})
else ()
  set(SOURCES ${SOURCES_COMPILER_VEC})
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
