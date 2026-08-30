set(SOURCES
    dip.f90
    em_field.f90
    initialization_rt.f90
    main_tddft.f90
    md.f90
    dg/rt_dg_overlapping_wannier.f90
    dg/rt_dg_hybrid_sparse_exchange.f90
    dg/rt_dg_hybrid_metric_solver.f90
    dg/rt_dg_hybrid_length_gauge.f90
    dg/rt_dg_hybrid_checkpoint.f90
    dg/rt_dg_hybrid_initialization.f90
    dg/rt_dg_hybrid_density_update.f90
    rt_angular_momentum.f90
    rt_local_chern_marker.f90
    rt_local_chern_marker_soi.f90
    plusu_fragment_support.f90
    taylor.f90
    time_evolution_step.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
