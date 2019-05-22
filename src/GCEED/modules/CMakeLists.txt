set(SOURCES
    allocate_mat.f90
    allocate_psl.f90
    calc_density.f90
    change_order.f90
    copyVlocal.f90
    copy_psi_mesh.f90
    deallocate_mat.f90
    gradient2.f90
    gradient.f90
    init_sendrecv.f90
    laplacian2.f90
    new_world.f90
    readbox_rt.f90
    read_pslfile.f90
    scf_data.f90
    sendrecv.f90
    sendrecv_self.f90
    sendrecv_tmp.f90
    sendrecv_groupob.f90
    sendrecv_groupob_ngp.f90
    sendrecv_groupob_tmp.f90
    sendrecvh.f90
    writebox_rt.f90
    persistent_comm.f90
    structure_opt.f90
    )
  #  share_mesh_1d_old.f90

list_prepend(SOURCES ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND SALMON_SOURCES ${SOURCES})
set(SALMON_SOURCES ${SALMON_SOURCES} PARENT_SCOPE)
