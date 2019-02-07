set(SOURCES
    eigen_subdiag_${FILE_MATHLIB}.f90
    eigen_subdiag_periodic_${FILE_MATHLIB}.f90
    gscg.f90
    gscg_periodic.f90
    rmmdiis_eigen.f90
    rmmdiis_eigen_lapack.f90
    rmmdiis.f90
    subspace_diag.f90
    subspace_diag_periodic.f90
    ybcg.f90
    ybcg_periodic.f90
   )

add_library(${SALMON_GS_LIB} STATIC ${SOURCES})
#add_dependencies(${SALMON_GS_LIB} ${SALMON_COMMON_LIB} ${SALMON_IO_LIB} ${SALMON_PARALLEL_LIB} ${SALMON_POISSON_LIB} ${SALMON_XC_LIB})
add_dependencies(${SALMON_GS_LIB} ${SALMON_IO_LIB} ${SALMON_PARALLEL_LIB} ${SALMON_COMMON_LIB})
list(APPEND SALMON_LINK_LIBRARIES ${SALMON_GS_LIB})
set(SALMON_LINK_LIBRARIES ${SALMON_LINK_LIBRARIES} PARENT_SCOPE)
