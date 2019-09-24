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
   zstencil_core/typical/omp.f90
   zstencil_core/typical/seq.f90
   )

if (USE_OPT_EXPLICIT_VECTORIZATION)
  set(SOURCES
      ${SOURCES}
      zstencil_core/${SIMD_SET}/omp.c
      zstencil_core/${SIMD_SET}/seq.c
      )
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
