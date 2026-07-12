!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
!=======================================================================
! End-to-end preparation of a stationary SALMON nodal DG initial state.
!=======================================================================
module rt_dg_nodal_salmon_prepare
  use structures, only: s_dft_system, s_parallel_info, s_orbital, s_scalar, s_stencil, s_pp_grid
  use rt_dg_fragment_types, only: s_dg_fragment_rt
  use rt_dg_nodal_types, only: s_dg_nodal_state
  use rt_dg_nodal_salmon_adapter, only: initialize_nodal_from_dg_coefficients, &
                                        build_nodal_local_potential, get_nodal_stencil_coefficients
  use rt_dg_nodal_salmon_ground_state, only: relax_nodal_salmon_ground_state_mpi
  use communication, only: comm_is_root
  implicit none
  private
  public :: prepare_nodal_salmon_ground_state

contains

  subroutine prepare_nodal_salmon_ground_state(dg_frag,system,info,spsi,Vh,Vxc,Vpsl,stencil,ppg, &
                                                halo_width,relaxation_step,max_iteration,tolerance,state, &
                                                vlocal,kinetic_center,kinetic_offset,gradient_offset, &
                                                eigenvalues,max_residual,niteration)
    type(s_dg_fragment_rt), intent(in) :: dg_frag
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_orbital), intent(in) :: spsi
    type(s_scalar), intent(in) :: Vh,Vxc(:),Vpsl
    type(s_stencil), intent(in) :: stencil
    type(s_pp_grid), intent(in) :: ppg
    integer, intent(in) :: halo_width,max_iteration
    real(8), intent(in) :: relaxation_step,tolerance
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), allocatable, intent(out) :: vlocal(:,:,:,:),eigenvalues(:,:)
    real(8), intent(out) :: kinetic_center,kinetic_offset(4,3),gradient_offset(4,3)
    real(8), intent(out) :: max_residual
    integer, intent(out) :: niteration
    integer :: ik
    real(8) :: t0,t1

    call cpu_time(t0)
    call initialize_nodal_from_dg_coefficients(dg_frag,system,info,halo_width,state)
    call cpu_time(t1)
    if (comm_is_root(dg_frag%id)) write(*,'(1x,a,1pe13.5)') '[DG-NODAL-PREP] seed_seconds=',t1-t0
    call build_nodal_local_potential(dg_frag,state,Vh,Vxc,Vpsl,vlocal)
    call cpu_time(t0)
    call get_nodal_stencil_coefficients(stencil,kinetic_center,kinetic_offset,gradient_offset)
    if (.not. allocated(ppg%zekr_uV)) stop 'nodal DG: nonlocal velocity-gauge projectors are unavailable'
    ik=lbound(ppg%zekr_uV,3)
    allocate(eigenvalues(state%nstate,state%nspin))
    if (comm_is_root(dg_frag%id)) then
      write(*,'(1x,a,i0,a,i0)') '[DG-NODAL-PREP] entering complete-H GS: nstate=',state%nstate, &
        ' ngrid_core=',product(state%core_size)
      flush(6)
    end if
    call relax_nodal_salmon_ground_state_mpi(state,dg_frag,ppg,ik,vlocal,kinetic_center, &
                                              kinetic_offset,gradient_offset,relaxation_step, &
                                              max_iteration,tolerance,dg_frag%icomm,eigenvalues, &
                                              max_residual,niteration)
    call cpu_time(t1)
    if (comm_is_root(dg_frag%id)) write(*,'(1x,a,1pe13.5)') '[DG-NODAL-PREP] gs_seconds=',t1-t0
    if (.not. state%dg_ground_state_ready) stop 'nodal DG: complete-H ground-state preparation failed'
  end subroutine prepare_nodal_salmon_ground_state

end module rt_dg_nodal_salmon_prepare
