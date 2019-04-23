set(SOURCES
    init.f90
    init_wf.f90
    fd_coef.f90
    prep_ps.f90
    pp_postprocess.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
