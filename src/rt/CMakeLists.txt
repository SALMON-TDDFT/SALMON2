set(SOURCES
    tdcdft_lrc.f90
    tdcdft_elf.f90
    dip.f90
    em_field.f90
    initialization_rt.f90
    main_tddft.f90
    md.f90
    taylor.f90
    time_evolution_step.f90
   )

if(USE_HSE)
  list(APPEND SOURCES hse_ptcn_core.f90 hse_ptcn.f90)
endif()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
