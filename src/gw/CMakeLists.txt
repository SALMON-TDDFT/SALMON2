set(SOURCES
    gw_coulomb.f90
    gw_mtxel.f90
    gw_epsilon.f90
    gw_symmetry.f90
    gw_sigma_x.f90
    gw_sigma_cohsex.f90
    gw_sigma_gpp.f90
    gw_qp.f90
    gw_qp_output.f90
    gw_main.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
