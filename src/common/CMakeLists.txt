set(SOURCES
   broyden.f90
   calc_allob.f90
   calc_iquotient.f90
   calc_iroot.f90
   calc_myob.f90
   calc_ob_num.f90
   check_corrkob.f90
   hpsi.f90
   inner_product.f90
   inner_product3.f90
   inner_product4.f90
   inner_product5.f90
   inner_product7.f90
   set_isstaend.f90
   stencil.f90
   update_overlap.f90
   structures.f90
   density_matrix.f90
   pack_unpack.f90
   sendrecv_grid.f90
   )

add_library(${SALMON_COMMON_LIB} STATIC ${SOURCES})
add_dependencies(${SALMON_COMMON_LIB} ${SALMON_MATH_LIB} ${SALMON_IO_LIB} ${SALMON_PARALLEL_LIB} ${SALMON_MISC_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_COMMON_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
