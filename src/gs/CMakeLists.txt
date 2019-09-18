set(SOURCES
    eigen_subdiag_${FILE_MATHLIB}.f90
    conjugate_gradient.f90
    mixing.f90
    rmmdiis.f90
    subspace_diagonalization.f90
    gram_schmidt_orth.f90
    scf_iteration.f90
    structure_opt.f90
    read_gs.f90
    band.f90
    occupation.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
