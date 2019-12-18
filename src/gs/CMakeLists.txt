set(SOURCES
    band.f90
    band_dft.f90
    conjugate_gradient.f90
    eigen_subdiag_${FILE_MATHLIB}.f90
    gram_schmidt_orth.f90
    init_gs.f90
    initialization_dft.f90
    main_dft.f90
    main_dft_md.f90
    mixing.f90
    scf_iteration_dft.f90
    scf_iteration.f90
    structure_opt.f90
    subspace_diagonalization.f90
    read_gs.f90
    occupation.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
