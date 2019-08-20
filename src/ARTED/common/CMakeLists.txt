set(SOURCES
    total_energy.f90
    ion_force.f90
    Exc_Cor.f90
    Hartree.f90
    psi_rho.f90
    hpsi${HPSI_TEST_SUFFIX}.f90
    restart.f90
    density_plot.f90
    projector.f90
    Ac_local_field.f90
    )

if (ENABLE_OPENACC)
  set(SOURCES ${SOURCES} acc/hpsi.f90)
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
