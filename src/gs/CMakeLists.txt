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
    gram_schmidt_orth.f90
    scf_iteration.f90
    structure_opt.f90
    read_gs.f90
    band.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
