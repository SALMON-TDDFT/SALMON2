set(SOURCES
    gram_schmidt_so.f90
    nondiagonal_so.f90
    prep_pp_so.f90
    pseudo_pt_so.f90
    spin_orbit_global.f90
    subspace_diagonalization_so.f90
    update_kvector_so.f90
    occupation_so.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
