set(SOURCES
    main_dft.f90
    main_dft_md.f90
    initialization_dft.f90
    scf_iteration_dft.f90
    band_dft.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
