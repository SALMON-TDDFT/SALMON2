set(SOURCES
    builtin_pz.f90
    builtin_pzm.f90
    builtin_pz_sp.f90
    builtin_pw.f90
    builtin_pw92.f90
    builtin_r2scan.f90
    builtin_tbmbj.f90
    salmon_xc.f90
   )

if(USE_HSE)
  list(APPEND SOURCES hse_wannier_gauge.f90 hse_wannier.f90 hse_symmetry.f90 hse_exchange.f90 hse_ace.f90 hse_semilocal.f90 hse_native.f90 hse_lcfo_rt.f90 lcfo_rt_wannier.f90 lcfo_ace_local.f90)
endif()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
