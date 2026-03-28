module rt_dg_fragment_parallel
  implicit none

contains

  subroutine setup_fragment_parallel_grid(dg_frag, mg_parent, info_frag, mg_frag)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use init_communicator, only: init_communicator_dft
    use initialization_sub, only: init_parallel_dft, init_grid_parallel
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_rgrid), intent(in) :: mg_parent
    type(s_parallel_info), intent(out) :: info_frag
    type(s_rgrid), intent(out) :: mg_frag

    type(s_dft_system) :: system_frag
    type(s_rgrid) :: lg_frag

    if (any(mg_parent%num(1:3) /= dg_frag%lgnum_total(1:3))) then
      stop 'RT-DG fragment-local MPI Hartree requires replicated parent real-space grid'
    end if

    system_frag%if_real_orbital = .true.
    system_frag%ngrid = product(mg_parent%num(1:3))
    system_frag%nspin = 1
    system_frag%no = 1
    system_frag%nk = 1
    system_frag%nion = 0
    system_frag%hvol = 0.0d0
    system_frag%hgs = 0.0d0
    system_frag%primitive_a = 0.0d0
    system_frag%det_a = 0.0d0
    system_frag%primitive_b = 0.0d0
    system_frag%rmatrix_a = 0.0d0
    system_frag%rmatrix_b = 0.0d0
    system_frag%mu = 0.0d0
    system_frag%vec_Ac = 0.0d0
    system_frag%vec_Ac_ext = 0.0d0
    system_frag%vec_E = 0.0d0
    system_frag%vec_E_ext = 0.0d0
    allocate(system_frag%vec_k(3, 1))
    system_frag%vec_k(:,1) = 0.0d0
    allocate(system_frag%wtk(1))
    system_frag%wtk(1) = 1.0d0
    allocate(system_frag%rocc(1, 1, 1))
    system_frag%rocc = 0.0d0

    info_frag%npk = 1
    info_frag%nporbital = 1
    info_frag%nprgrid(1:3) = nproc_rgrid(1:3)
    call init_communicator_dft(dg_frag%icomm_frag, info_frag)
    call init_parallel_dft(system_frag, info_frag)

    lg_frag = mg_parent
    call init_grid_parallel(info_frag, lg_frag, mg_frag)

    deallocate(system_frag%vec_k)
    deallocate(system_frag%wtk)
    deallocate(system_frag%rocc)
  end subroutine setup_fragment_parallel_grid

  subroutine setup_fragment_system(dg_frag, system, info, mg_parent, system_frag, info_frag, mg_frag)
    use structures
    use rt_dg_fragment_types, only: s_dg_fragment_rt
    use init_communicator, only: init_communicator_dft
    use initialization_sub, only: init_parallel_dft, init_grid_parallel
    use salmon_global, only: nproc_rgrid
    implicit none
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_rgrid), intent(in) :: mg_parent
    type(s_dft_system), intent(out) :: system_frag
    type(s_parallel_info), intent(out) :: info_frag
    type(s_rgrid), intent(out) :: mg_frag

    integer :: no_src, no_basis_global
    type(s_rgrid) :: lg_frag

    if (info%npk /= 1 .or. info%nporbital /= 1 .or. product(info%nprgrid) /= 1) then
      stop 'RT-DG fragment-local MPI basis update stage-1 requires parent k/orbital/r-space replication'
    end if
    if (any(mg_parent%num(1:3) /= dg_frag%lgnum_total(1:3))) then
      stop 'RT-DG fragment-local MPI basis update stage-1 requires replicated parent real-space grid'
    end if

    system_frag = system
    no_basis_global = max(1, maxval(dg_frag%n_basis(1:dg_frag%n_frag,1:dg_frag%nspin)))
    system_frag%no = min(system%no, no_basis_global)
    system_frag%nk = 1
    system_frag%if_real_orbital = .true.

    if (allocated(system_frag%vec_k)) deallocate(system_frag%vec_k)
    allocate(system_frag%vec_k(3, system_frag%nk))
    system_frag%vec_k(:,1) = 0.0d0

    if (allocated(system_frag%wtk)) deallocate(system_frag%wtk)
    allocate(system_frag%wtk(system_frag%nk))
    system_frag%wtk(1) = 1.0d0

    no_src = 0
    if (allocated(system%rocc)) no_src = size(system%rocc,1)
    if (allocated(system_frag%rocc)) deallocate(system_frag%rocc)
    allocate(system_frag%rocc(system_frag%no, system_frag%nk, system_frag%nspin))
    system_frag%rocc = 0.0d0
    if (allocated(system%rocc) .and. no_src > 0) then
      system_frag%rocc(1:min(system_frag%no,no_src),1,1:system_frag%nspin) = &
        system%rocc(1:min(system_frag%no,no_src),1,1:system_frag%nspin)
    end if

    info_frag%npk = 1
    info_frag%nporbital = 1
    info_frag%nprgrid(1:3) = nproc_rgrid(1:3)
    call init_communicator_dft(dg_frag%icomm_frag, info_frag)
    call init_parallel_dft(system_frag, info_frag)

    lg_frag = mg_parent
    call init_grid_parallel(info_frag, lg_frag, mg_frag)

  end subroutine setup_fragment_system

  subroutine finalize_fragment_parallel(info_frag)
    use communication, only: comm_free_group, COMM_GROUP_NULL
    use structures, only: s_parallel_info
    implicit none
    type(s_parallel_info), intent(inout) :: info_frag

    if (allocated(info_frag%irank_io)) deallocate(info_frag%irank_io)
    if (allocated(info_frag%io_s_all)) deallocate(info_frag%io_s_all)
    if (allocated(info_frag%io_e_all)) deallocate(info_frag%io_e_all)
    if (allocated(info_frag%numo_all)) deallocate(info_frag%numo_all)
    if (allocated(info_frag%imap)) deallocate(info_frag%imap)
    if (allocated(info_frag%imap_isolated_ffte)) deallocate(info_frag%imap_isolated_ffte)

    if (info_frag%icomm_xy /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_xy)
    if (info_frag%icomm_z /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_z)
    if (info_frag%icomm_y /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_y)
    if (info_frag%icomm_x /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_x)
    if (info_frag%icomm_ko /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_ko)
    if (info_frag%icomm_ro /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_ro)
    if (info_frag%icomm_k /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_k)
    if (info_frag%icomm_o /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_o)
    if (info_frag%icomm_r /= COMM_GROUP_NULL) call comm_free_group(info_frag%icomm_r)
  end subroutine finalize_fragment_parallel

end module rt_dg_fragment_parallel
