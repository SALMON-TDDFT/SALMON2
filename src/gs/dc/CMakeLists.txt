set(SOURCES
    dc_fragment_geometry.f90
    dcdft.f90
    dcdft_soi.f90
    lcfo.f90
    lcfo_soi.f90
    lcfo_soi_init.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
