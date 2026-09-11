set(SOURCES
    dip.f90
    em_field.f90
    initialization_rt.f90
    main_tddft.f90
    md.f90
    rt_angular_momentum.f90
    rt_local_chern_marker.f90
    rt_local_chern_marker_soi.f90
    plusu_fragment_support.f90
    taylor.f90
    time_evolution_step.f90
   )

if (USE_MPI AND USE_SCALAPACK)
  list(APPEND SOURCES
    dg/rt_dg_hybrid_sparse_exchange.f90
    dg/rt_dg_hybrid_point_density.f90
    dg/rt_dg_hybrid_metric_solver.f90
    dg/rt_dg_hybrid_length_gauge.f90
    dg/rt_dg_hybrid_system_identity.f90
    dg/rt_dg_hybrid_checkpoint_v5.f90
    dg/rt_dg_hybrid_refinement_receipt.f90
    dg/rt_dg_hybrid_occupied_checkpoint.f90
    dg/rt_dg_hybrid_structural_graph.f90
    dg/rt_dg_hybrid_sparse_projection.f90
    dg/rt_dg_hybrid_initialization_v5.f90
    dg/rt_dg_hybrid_density_update.f90
    dg/rt_dg_hybrid_stationarity.f90
   )
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
