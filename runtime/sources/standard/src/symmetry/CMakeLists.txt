set(SOURCES
    symmetry.f90
    symmetry_kvector.f90
    symmetry_rho.f90
    symmetry_vector.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
