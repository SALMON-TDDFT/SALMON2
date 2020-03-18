set(SOURCES
   broyden.f90
   const.f90
   hamiltonian.f90
   initialization.f90
   stencil.f90
   nonlocal_potential.f90
   structures.f90
   density_matrix.f90
   pack_unpack.f90
   sendrecv_grid.f90
   total_energy.f90
   lattice.f90
   force.f90
   calc_gradient_psi.f90
   )

 #if (IS_FUJITSU_COMPILER)
  set(SOURCES ${SOURCES}
      zstencil_core/typical/seq_frt.f90
      zstencil_core/typical/omp_frt.f90
      zstencil_microac_core/typical/seq_frt.f90
      zstencil_microac_core/typical/omp_frt.f90
      )
    #else ()
    #  set(SOURCES ${SOURCES}
    #      zstencil_core/typical/seq.f90
    #      zstencil_core/typical/omp.f90
    #      zstencil_microac_core/typical/seq.f90
    #      zstencil_microac_core/typical/omp.f90
    #      )
    #endif ()

if (USE_OPT_EXPLICIT_VECTORIZATION)
  set(SOURCES
      ${SOURCES}
      zstencil_core/${SIMD_SET}/omp.c
      zstencil_core/${SIMD_SET}/seq.c
      )
endif ()

if (USE_SCALAPACK)
  set(SOURCES ${SOURCES}
      scalapack_module.f90
      )
endif ()

if (USE_EIGENEXA)
  set(SOURCES ${SOURCES}
      eigenexa_module.f90
      )
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
