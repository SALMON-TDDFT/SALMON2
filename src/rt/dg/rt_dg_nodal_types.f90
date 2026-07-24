!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
! Thin compatibility layer for the historical RT nodal API.
module rt_dg_nodal_types
  use dg_nodal_state, only: s_dg_nodal_state => s_dg_nodal_common_state
  use dg_nodal_state, only: nodal_face_info => s_dg_nodal_face
  use dg_nodal_state, only: nodal_face_count
  use dg_nodal_state, only: dg_nodal_face_slot
  use dg_nodal_state, only: initialize_dg_nodal_common_state_local
  use dg_nodal_state, only: allocate_dg_nodal_face_buffers
  implicit none
  private
  public :: s_dg_nodal_state, nodal_face_info, nodal_face_count
  public :: nodal_face_slot, allocate_nodal_state, allocate_nodal_face_buffers
  public :: accept_nodal_dg_ground_state

contains

  pure integer function nodal_face_slot(axis,side) result(iface)
    integer, intent(in) :: axis,side
    iface=dg_nodal_face_slot(axis,side)
  end function nodal_face_slot

  subroutine allocate_nodal_state(state,fragment,core_size,halo_width,nstate,nspin)
    type(s_dg_nodal_state), intent(inout) :: state
    integer, intent(in) :: fragment,core_size(3),halo_width,nstate,nspin
    integer(8), allocatable :: owner(:,:,:)
    logical :: ok
    character(len=256) :: message
    integer :: i

    allocate(owner(core_size(1),core_size(2),core_size(3)))
    owner=reshape([(int(i,8),i=1,size(owner))],shape(owner))
    call initialize_dg_nodal_common_state_local(state,fragment,core_size,core_size,[0,0,0],halo_width, &
      nstate,nspin,owner,int(fragment,8),int(fragment,8),ok,message)
    deallocate(owner)
    if (.not.ok) error stop trim(message)
  end subroutine allocate_nodal_state

  subroutine allocate_nodal_face_buffers(face,core_size,halo_width,nstate,nspin)
    type(nodal_face_info), intent(inout) :: face
    integer, intent(in) :: core_size(3),halo_width,nstate,nspin
    logical :: ok
    call allocate_dg_nodal_face_buffers(face,core_size,halo_width,nstate,nspin,ok)
    if (.not.ok) error stop 'nodal DG: invalid face buffer request'
  end subroutine allocate_nodal_face_buffers

  subroutine accept_nodal_dg_ground_state(state,residual,tolerance)
    type(s_dg_nodal_state), intent(inout) :: state
    real(8), intent(in) :: residual,tolerance
    if (residual < 0.0d0 .or. residual /= residual) error stop 'nodal DG: invalid ground-state residual'
    if (tolerance <= 0.0d0) error stop 'nodal DG: invalid ground-state tolerance'
    state%dg_ground_state_residual=residual
    if (residual > tolerance) error stop 'nodal DG: ground state is not stationary under the DG Hamiltonian'
    state%ground_state_ready=.true.
  end subroutine accept_nodal_dg_ground_state
end module rt_dg_nodal_types
