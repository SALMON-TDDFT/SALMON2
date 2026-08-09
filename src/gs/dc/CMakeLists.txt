set(SOURCES
    dcdft.f90
    lcfo.f90
   )

if (USE_SCALAPACK)
  list(INSERT SOURCES 1 lcfo_diag_chefsi.f90)
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
