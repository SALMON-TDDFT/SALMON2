set(SOURCES
    dc_fragment_geometry.f90
    dc_scf_convergence.f90
    dc_fragment_occupation.f90
    dg_dc_seed_checkpoint.f90
    dg_overlapping_wannier_types.f90
    dg_spatial_mlwf_seed.f90
    dg_hybrid_wannier_selection.f90
    dg_hybrid_full_cell_operator_adapter.f90
    dg_hybrid_occupation_policy.f90
    dg_hybrid_ground_state_types.f90
    dg_hybrid_continuation_state.f90
    dg_hybrid_continuation_controller.f90
    dg_hybrid_low_energy_symmetry.f90
    dg_hybrid_sipg_operator.f90
    dg_hybrid_broken_volume.f90
    dg_hybrid_variational_payload.f90
    dg_hybrid_divided_operator.f90
    dg_hybrid_real_space_residual.f90
    dg_hybrid_variational_potential.f90
    dg_hybrid_generalized_eigensystem.f90
    dg_hybrid_density.f90
    dg_hybrid_scf.f90
    dg_hybrid_fragment_basis.f90
    dg_hybrid_production_face_traces.f90
    dg_hybrid_fragment_basis_stream.f90
    dg_hybrid_projected_fragment_pipeline.f90
    dg_hybrid_production_fragment_basis.f90
    dg_hybrid_fragment_solver.f90
    dg_hybrid_divided_scf.f90
    dg_hybrid_lcfo_support_redistribution.f90
    dg_hybrid_lcfo.f90
    dg_hybrid_block_cg.f90
    dg_hybrid_fragment_subspace.f90
    dg_hybrid_fragment_preconditioner.f90
    dg_overlapping_wannier_symmetry.f90
    dg_overlapping_wannier_localization.f90
    dg_overlapping_wannier_projection.f90
    dg_overlapping_wannier_metric.f90
    dg_overlapping_wannier_w90.f90
    dg_hybrid_fragment_wannier.f90
    dg_hybrid_fragment_selection.f90
    dg_hybrid_fragment_admission.f90
    dg_hybrid_production_support.f90
    dg_hybrid_fragment_thermal.f90
    dg_overlapping_wannier_construction.f90
    dg_hybrid_localization_first.f90
    dg_hybrid_certified_rt_basis.f90
    dg_hybrid_retained_basis_symmetry.f90
    dg_overlapping_wannier_operators.f90
    dg_overlapping_wannier_nonlocal.f90
    dg_nonlocal_projector_range.f90
    dg_overlapping_wannier_observables.f90
    dg_overlapping_wannier_full_cell.f90
    dg_overlapping_wannier_solver.f90
    dg_overlapping_wannier_density.f90
    dg_overlapping_wannier_scf.f90
    dg_overlapping_wannier_checkpoint.f90
  dg_dc_buffer_core_faces.f90
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
