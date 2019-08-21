set(SOURCES
    calc_dos.f90
    calc_occupation.f90
    calc_pdos.f90
    check_dos_pdos.f90
    convert_input_scf.f90
    copy_density.f90
    deallocate_sendrecv.f90
    prep_ini.f90
    real_space_dft.f90
    subdgemm_lapack.f90
    write_eigen.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
