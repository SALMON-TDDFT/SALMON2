set(SOURCES
    hartree.f90
    init_poisson.f90
    poisson_dirichlet.f90
    poisson_isolated.f90
    poisson_periodic.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
