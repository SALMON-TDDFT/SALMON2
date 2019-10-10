set(SOURCES
    calc_elf.f90
    set_filename.f90
    set_icoo1d.f90
    write_psi.f90
    ylm.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
