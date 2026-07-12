set(SOURCES
    dc_fragment_geometry.f90
    dcdft.f90
    dcdft_soi.f90
    lcfo_wannier_sawf.f90
    lcfo_wannier_sawf_basis.f90
    lcfo_wannier_sawf_band.f90
    lcfo_wannier_sawf_collective.f90
    lcfo_wannier_sawf_dmn.f90
    lcfo_wannier_sawf_win.f90
    lcfo_wannier_command.f90
    lcfo.f90
    lcfo_flux.f90
    lcfo_soi.f90
   lcfo_soi_init.f90
   )

if (HAVE_SPGLIB)
  list(APPEND SOURCES lcfo_wannier_spglib.c)
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
