set(SOURCES
    hartree.f90
    hartree_boundary.f90
    hartree_cg.f90
    hartree_ffte.f90
    hartree_periodic.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
