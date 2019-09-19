set(SOURCES
    allocate_mat.f90
    deallocate_mat.f90
    gradient.f90
    new_world.f90
    scf_data.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
