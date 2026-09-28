set(SOURCES
    dc_thermal.f90
    dcdft.f90
    dc_force.f90
    dc_projector_force.f90
    lcfo_mesh_stream.f90
    lcfo_mesh_tile.f90
    lcfo_complex.f90
    lcfo.f90
   )

if (USE_SCALAPACK)
  list(APPEND SOURCES lcfo_scalapack.f90)
  list(INSERT SOURCES 1 lcfo_diag_chefsi.f90)
  list(INSERT SOURCES 2 lcfo_diag_chefsi_complex.f90)
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
