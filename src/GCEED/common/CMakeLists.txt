set(SOURCES
    allocate_sendrecv.f90
    calcELF.f90
    calcJxyz_all.f90
    calcJxyz_all_periodic.f90
    calc_pmax.f90
    calcuV.f90
    calcVpsl.f90
    calcVpsl_periodic.f90
    calcVpsl_periodic_FFTE.f90
    check_cep.f90
    check_fourier.f90
    conv_p0.f90
    conv_p.f90
    exc_cor_ns.f90
    init_wf.f90
    init_k_rd.f90
    OUT_IN_data.f90
    psl.f90
    prep_poisson_fft.f90
    read_copy_pot.f90
    set_filename.f90
    set_icoo1d.f90
    set_imesh_oddeven.f90
    set_is.f90
    storevpp.f90
    writepsi.f90
    ylm.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
