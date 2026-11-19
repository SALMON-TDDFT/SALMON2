set(SOURCES
    calc_jxyz_plusu.f90
    density_matrix_and_energy_plusu.f90
    prep_pp_plusu.f90
    update_kvector_plusu.f90
    calc_nps_plusu.f90
    plusu_global.f90
    pseudo_pt_plusu.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
