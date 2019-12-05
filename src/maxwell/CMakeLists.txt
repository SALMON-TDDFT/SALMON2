set(SOURCES
    salmon_maxwell.f90
    main_maxwell.f90
    weyl_init.f90
    weyl_finalize.f90
    weyl_calc.f90
    fdtd_coulomb_gauge.f90
    set_coo_em.f90
    find_point_em.f90
    eh_calc.f90
    eh_init.f90
    eh_finalize.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
