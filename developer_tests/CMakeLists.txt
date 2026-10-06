# Optional CPU numerical development checks; separate from input-based regressions.
if(NOT USE_OPENACC)
  add_subdirectory(652_dc_lcfo/rt)
  if(USE_HSE)
    add_subdirectory(651_hybrid_exchange/history)
    add_subdirectory(651_hybrid_exchange/lapack)
  endif()
endif()
