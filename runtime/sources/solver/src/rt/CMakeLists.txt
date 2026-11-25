set(SOURCES
    tddft_worker.f90
    dip.f90
    em_field.f90
    near_field_coupling.f90
    grid_field.f90
    orbital_response.f90
    solid_harmonics.f90
    initialization_rt.f90
    main_tddft.f90
    md.f90
    taylor.f90
    step_recovery.f90
    time_evolution_step.f90
   )

if(USE_HSE)
  list(APPEND SOURCES exx_ptcn.f90)
endif()

if(USE_HSE OR USE_SCALAPACK)
  list(APPEND SOURCES ptcn_core.f90 exx_ptcn_core.f90 ptcn_native.f90
    ptcn_rectangular_fft.f90 ptcn_kinetic_preconditioner.f90)
endif()
if(USE_SCALAPACK)
  list(APPEND SOURCES ptcn_distributed_matrix.f90 ptcn_native_matrix.f90)
endif()

if(USE_MPI)
  list(APPEND SOURCES embedded_tddft.f90 ../salmonpy/parallel_bridge.f90)
endif()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
