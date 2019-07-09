set(SOURCES
    add_polynomial.f90
    calc_current.f90
    calc_env_trigon.f90
    calcAext.f90
    calcEstatic.f90
    calcVbox.f90
    check_ae_shape.f90
    convert_input_rt.f90
    dip.f90
    gradient_ex.f90
    initA.f90
    projection.f90
    read_rt.f90
    real_time_dft.f90
    set_vonf_sd.f90
    taylor_coe.f90
    time_evolution_step.f90
    xc_fast.f90
    )
  #    aft_fourier.f90 # Is this single program?

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
