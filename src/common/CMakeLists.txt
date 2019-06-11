set(SOURCES
   broyden.f90
   calc_allob.f90
   calc_iquotient.f90
   calc_iroot.f90
   calc_myob.f90
   calc_ob_num.f90
   check_corrkob.f90
   check_mg.f90
   check_ng.f90
   check_numcpu.f90
   const.f90
   hpsi.f90
   inner_product.f90
   inner_product3.f90
   inner_product4.f90
   inner_product5.f90
   inner_product7.f90
   set_isstaend.f90
   stencil.f90
   pseudo_pt.f90
   update_overlap.f90
   set_numcpu.f90
   structures.f90
   density_matrix.f90
   pack_unpack.f90
   sendrecv_grid.f90
   total_energy.f90
   writefield.f90
   writefile3d.f90
   lattice.f90
   force.f90
   print.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
