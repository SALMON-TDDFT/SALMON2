set(SOURCES
    calc_allob.f90
    calc_iquotient.f90
    calc_iroot.f90
    calc_myob.f90
    calc_ob_num.f90
    check_corrkob.f90
    calc_elf.f90
    calc_pmax.f90
    calcVpsl_periodic_FFTE.f90
    check_fourier.f90
    conv_p0.f90
    conv_p.f90
    init_wf.f90
    init_k_rd.f90
    read_write_gs_bin.f90
    read_copy_pot.f90
    set_filename.f90
    set_icoo1d.f90
    set_imesh_oddeven.f90
    write_psi.f90
    ylm.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
