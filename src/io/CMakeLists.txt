set(SOURCES
    inputoutput.f90
    salmon_file.f90
    salmon_global.f90
    flops.f90
    perflog.f90
   )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
