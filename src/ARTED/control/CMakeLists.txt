set(SOURCES
    control_sc.f90
    control_ms.f90
    control_ms_raman.f90
    inputfile.f90
    initialization.f90
    ground_state.f90
    optimization.f90
    md_gs.f90
    ana_RT_useGS.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
