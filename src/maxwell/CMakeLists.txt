set(SOURCES
    common_maxwell.f90
    fdtd_eh.f90
    main_maxwell.f90
    fdtd_coulomb_gauge.f90
    fdtd_weyl_gauge.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
