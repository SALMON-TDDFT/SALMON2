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
set(SOURCES_ACC
    F90/acc/current.f90
    F90/acc/hpsi.f90
    F90/acc/total_energy.f90
    )

###
if (EXPLICIT_VEC)
  set(SOURCES_OMP ${SOURCES_EXPLICIT_VEC})
else ()
  set(SOURCES_OMP ${SOURCES_COMPILER_VEC})
endif ()

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES_ACC} ${SOURCES_COMPILER_VEC})
else ()
  set(SOURCES ${SOURCES_OMP})
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
