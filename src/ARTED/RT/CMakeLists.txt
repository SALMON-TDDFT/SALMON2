set(SOURCES
    dt_evolve.f90
    current.f90
    init_Ac.f90
    Fourier_tr.f90
    hamiltonian.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
