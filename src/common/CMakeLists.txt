set(SOURCES
   broyden.f90
   calc_allob.f90
   calc_iquotient.f90
   calc_iroot.f90
   calc_myob.f90
   calc_ob_num.f90
   check_corrkob.f90
   check_numcpu.f90
   const.f90
   hpsi.f90
   initialization.f90
   inner_product.f90
   inner_product3.f90
   set_isstaend.f90
   stencil.f90
   pseudo_pt.f90
   set_numcpu.f90
   structures.f90
   density_matrix.f90
   pack_unpack.f90
   sendrecv_grid.f90
   set_gridcoordinate.f90
   total_energy.f90
   write_field.f90
   write_file3d.f90
   write.f90
   lattice.f90
   force.f90
   zstencil_core/typical/omp.f90
   zstencil_core/typical/seq.f90
   )

if (EXPLICIT_VEC)
  set(SOURCES
      ${SOURCES}
      zstencil_core/${SIMD_SET}/omp.c
      zstencil_core/${SIMD_SET}/seq.c
      )
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
