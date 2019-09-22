set(SOURCES
    taylor.f90
    read_write_restart_rt.f90
    md.f90
    em_field.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
