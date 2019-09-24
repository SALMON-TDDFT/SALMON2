set(SOURCES
    salmon_parallel.f90
    code_optimization.f90
    init_communicator.f90
    set_numcpu.f90
)

if (USE_MPI)
  set(SOURCES ${SOURCES} salmon_communication.f90)
else ()
  set(SOURCES ${SOURCES} salmon_communication_dummy.f90)
endif ()

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
