set(SOURCES
    builtin_pz.f90
    builtin_pzm.f90
    builtin_pz_sp.f90
    builtin_tbmbj.f90
    salmon_xc.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
