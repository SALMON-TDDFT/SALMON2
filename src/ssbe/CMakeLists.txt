set(SOURCES
    util_ssbe.f90
    mpiio_export_ssbe.f90
    datafile_ssbe.f90
    input_checker_sbe.f90
    gs_info_ssbe.f90
    common_ssbe.f90
    degenerate_block_ssbe.f90
    sbe_collision_gw_core.f90
    sbe_collision_gw.f90
    bloch_solver_ssbe.f90
    realtime_ssbe.f90
    shaper_ssbe.f90
    multiscale_ssbe.f90
    main_ssbe.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
