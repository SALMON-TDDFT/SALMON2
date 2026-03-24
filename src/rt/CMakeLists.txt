set(SOURCES
    dip.f90
    em_field.f90
    initialization_rt.f90
    main_tddft.f90
    md.f90
    dg/rt_dg_fragment_types.f90
    dg/rt_dg_basis_projection.f90
    dg/rt_dg_plane_wave.f90
    rt_angular_momentum.f90
    rt_local_chern_marker.f90
    rt_local_chern_marker_soi.f90
    dg/rt_dg_fragment_ops.f90
    dg/rt_dg_fragment.f90
    dg/rt_dg_fragment_soi.f90
	    plusu_fragment_support.f90
    taylor.f90
    time_evolution_step.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
