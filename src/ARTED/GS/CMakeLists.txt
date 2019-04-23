set(SOURCES
    io_gs_wfn_k.f90
    Gram_Schmidt.f90
    diag.f90
    sp_energy.f90
    Occupation_Redistribution.f90
    CG.f90
    write_GS_data.f90
    Fermi_Dirac_distribution.f90
    )

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
