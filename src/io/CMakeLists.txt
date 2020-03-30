set(SOURCES
    checkpoint_restart.f90
    inputoutput.f90
    salmon_global.f90
    perflog.f90
    posix.c
    filesystem.f90
    write_field.f90
    write_file3d.f90
    write.f90
    main_dft2tddft.f90
    main_dft_k_expand.f90
    main_dft_k_expand_slice.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
