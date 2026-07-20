set(SOURCES
    dc_fragment_geometry.f90
    dg_wpw_lda_hartree.f90
    dg_wannier_pw_scf.f90
    dg_wpw_matrix_free_scf.f90
    dg_wpw_scf_driver.f90
    dg_wpw_rank_local_quadrature.f90
    dg_wpw_core_point_builder.f90
    dg_wpw_core_w_provider.f90
    dg_wpw_occupied_w_basis.f90
    dg_wpw_lcfo_operator_adapter.f90
    dg_wpw_lcfo_potential_map.f90
    dg_wpw_lcfo_volume_operator.f90
    dg_wpw_wannier_tail_halo.f90
    dcdft.f90
    dcdft_soi.f90
    lcfo_wannier_sawf.f90
    lcfo_wannier_sawf_templates.f90
    lcfo_wannier_sawf_orchestrator.f90
    lcfo_wannier_sawf_seed.f90
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
