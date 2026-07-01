set(SOURCES
    dip.f90
    em_field.f90
    initialization_rt.f90
    main_tddft.f90
    md.f90
    dg/rt_dg_fragment_parallel.f90
    dg/rt_dg_fragment_types.f90
    dg/rt_dg_basis_projection.f90
    dg/rt_dg_wpw_linalg.f90
    dg/rt_dg_wpw_window.f90
    dg/rt_dg_wpw_reduced_mixedz_diag.f90
    dg/rt_dg_plane_wave.f90
    rt_angular_momentum.f90
    rt_local_chern_marker.f90
    rt_local_chern_marker_soi.f90
    dg/rt_dg_fragment_ops.f90
    dg/rt_dg_fragment_coefficients.f90
    dg/rt_dg_fragment_layout.f90
    dg/rt_dg_fragment_lifecycle.f90
    dg/rt_dg_fragment_io.f90
    dg/rt_dg_fragment_io_soi.f90
    dg/rt_dg_hse_exchange.f90
    dg/rt_dg_fragment.f90
    dg/rt_dg_fragment_hamiltonian_soi.f90
    dg/rt_dg_fragment_soi.f90
	    plusu_fragment_support.f90
    taylor.f90
    time_evolution_step.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)


set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment.f90
    PROPERTIES
    OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment_halo.f90;${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment_hamiltonian.f90;${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment_parallel.f90"
)

set_source_files_properties(
    ${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment_soi.f90
    PROPERTIES
    OBJECT_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment_halo_soi.f90;${CMAKE_CURRENT_SOURCE_DIR}/dg/rt_dg_fragment_parallel.f90"
)
