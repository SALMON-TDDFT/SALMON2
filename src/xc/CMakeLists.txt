set(SOURCES
    builtin_pz.f90
    builtin_pzm.f90
    builtin_pz_sp.f90
    builtin_pw.f90
    builtin_pw92.f90
    builtin_r2scan.f90
    builtin_tbmbj.f90
    exx_functional.f90
    exx_gs_metadata.f90
    salmon_xc.f90
   )

if(USE_HSE)
  list(APPEND SOURCES
    exx_k_backend.f90
    exx_k_cufft.f90
    exx_batch_backend.f90
    exx_cufft.f90
    exx_local_fft.f90
    exx_adaptive_support.f90
    exx_pair_candidates.f90
    exx_spatial_local.f90
    rvv10.f90
    fftw_pencils.f90
    rvv10_distributed.f90
    rvv10_native.f90
    rvv10_ffte.f
    exx_wannier_gauge.f90
    exx_wannier.f90
    exx_symmetry.f90
    exx_k_exchange.f90
    exx_ace.f90
    exx_distributed_metric.f90
    exx_orbitals.f90
    exx_distributed_gauge.f90
    hybrid_semilocal.f90
    exx_spatial.f90
    exx_native.f90
    exx_lcfo_rt.f90
    lcfo_wf_support.f90
    lcfo_projection.f90
    lcfo_rt_wannier.f90
    lcfo_ace_local.f90
    lcfo_gram.f90
    lcfo_mlwf_links.f90
    lcfo_seed.f90
    lcfo_dist_rows.f90
    lcfo_dist_dense.f90
  )
endif()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
