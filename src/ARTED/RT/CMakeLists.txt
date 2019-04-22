set(SOURCES
    dt_evolve.f90
    current.f90
    init_Ac.f90
    Fourier_tr.f90
    )

set(SOURCES_ACC
    acc/current.f90
    acc/hamiltonian.f90
    )

set(SOURCES_OMP
    hamiltonian.f90
    )

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES} ${SOURCES_ACC})
else ()
  set(SOURCES ${SOURCES} ${SOURCES_OMP})
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
