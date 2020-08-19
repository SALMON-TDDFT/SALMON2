set(SOURCES
    pp/read_ps_upf_module.f90
    pp/input_pp.f90
    pp/salmon_pp.f90
    pp/prep_pp.f90
    jellium.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
