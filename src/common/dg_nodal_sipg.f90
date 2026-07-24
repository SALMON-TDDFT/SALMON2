!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_nodal_sipg
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  type, public :: s_dg_nodal_sipg_face
    integer(8) :: global_face_id=0
    integer(8) :: trace_epoch=0
    integer :: fragment_minus=0, fragment_plus=0, axis=0
    integer :: periodic_shift(3)=0
    logical :: canonical_owner=.false.
    logical :: physical_face=.false.
    logical :: physical_supercell_periodic=.false.
    logical :: auxiliary_periodic_wrap=.false.
    real(8) :: normal(3)=0.0d0
    real(8) :: h_normal=0.0d0, face_weight=0.0d0
    complex(8) :: u_minus=(0.0d0,0.0d0), u_plus=(0.0d0,0.0d0)
    complex(8) :: dn_minus=(0.0d0,0.0d0), dn_plus=(0.0d0,0.0d0)
  end type s_dg_nodal_sipg_face

  type, public :: s_dg_nodal_sipg_action
    complex(8) :: consistency_value(2)=(0.0d0,0.0d0)
    complex(8) :: symmetry_normal(2)=(0.0d0,0.0d0)
    complex(8) :: penalty_value(2)=(0.0d0,0.0d0)
    complex(8) :: total_value(2)=(0.0d0,0.0d0)
    complex(8) :: total_normal(2)=(0.0d0,0.0d0)
  end type s_dg_nodal_sipg_action

  type, public :: s_dg_nodal_sipg_diagnostics
    real(8) :: hermiticity_defect=0.0d0
    real(8) :: internal_cancellation_defect=0.0d0
    real(8) :: jump_norm=0.0d0
    real(8) :: penalty_energy=0.0d0
    integer(8) :: trace_epoch=0
    integer :: ownership_count=0
  end type s_dg_nodal_sipg_diagnostics

  public :: evaluate_dg_nodal_sipg_face, apply_dg_nodal_sipg_faces_mpi

contains

  pure subroutine evaluate_dg_nodal_sipg_face(u_minus,u_plus,dn_minus,dn_plus,h_normal,face_weight, &
                                               penalty_factor,action,info)
    complex(8), intent(in) :: u_minus,u_plus,dn_minus,dn_plus
    real(8), intent(in) :: h_normal,face_weight,penalty_factor
    type(s_dg_nodal_sipg_action), intent(out) :: action
    integer, intent(out) :: info
    complex(8) :: jump,average_dn
    real(8) :: penalty_over_h

    action=s_dg_nodal_sipg_action()
    info=1
    if (h_normal <= 0.0d0 .or. face_weight < 0.0d0 .or. penalty_factor <= 0.0d0) return
    if (.not.ieee_is_finite(h_normal) .or. .not.ieee_is_finite(face_weight) .or. &
        .not.ieee_is_finite(penalty_factor)) return
    if (.not.finite_complex(u_minus) .or. .not.finite_complex(u_plus) .or. &
        .not.finite_complex(dn_minus) .or. .not.finite_complex(dn_plus)) return
    jump=u_minus-u_plus
    average_dn=0.5d0*(dn_minus+dn_plus)
    penalty_over_h=penalty_factor/h_normal
    action%consistency_value=face_weight*[-average_dn,average_dn]
    action%symmetry_normal=face_weight*[-0.5d0*jump,-0.5d0*jump]
    action%penalty_value=face_weight*[penalty_over_h*jump,-penalty_over_h*jump]
    action%total_value=action%consistency_value+action%penalty_value
    action%total_normal=action%symmetry_normal
    info=0
  end subroutine evaluate_dg_nodal_sipg_face

  subroutine apply_dg_nodal_sipg_faces_mpi(faces,expected_epoch,penalty_factor,communicator, &
                                            action,diagnostics,ok,message)
    type(s_dg_nodal_sipg_face), intent(in) :: faces(:)
    integer(8), intent(in) :: expected_epoch
    real(8), intent(in) :: penalty_factor
    integer, intent(in) :: communicator
    type(s_dg_nodal_sipg_action), allocatable, intent(inout) :: action(:)
    type(s_dg_nodal_sipg_diagnostics), intent(inout) :: diagnostics
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_sipg_action), allocatable :: candidate(:)
    type(s_dg_nodal_sipg_diagnostics) :: candidate_diagnostics
    integer :: i,info
    logical :: local_ok
    real(8) :: local_sums(3),global_sums(3),local_energy_imag,local_balance
    complex(8) :: energy
    type(s_dg_nodal_sipg_action) :: reverse_action
#ifdef USE_MPI
    integer :: ierr
    logical :: stage_ok
#endif

    allocate(candidate(size(faces)))
    candidate=s_dg_nodal_sipg_action()
    local_ok=expected_epoch > 0_8 .and. penalty_factor > 0.0d0 .and. ieee_is_finite(penalty_factor)
    do i=1,size(faces)
      local_ok=local_ok .and. valid_face(faces(i),expected_epoch)
    end do
    call collective_and(local_ok,communicator,ok)
    if (.not.ok) then
      message='nodal SIPG: collective face validation failed'
      deallocate(candidate)
      return
    end if

    do i=1,size(faces)
      if (faces(i)%canonical_owner .and. faces(i)%physical_face .and. &
          .not.faces(i)%auxiliary_periodic_wrap) then
        call evaluate_dg_nodal_sipg_face(faces(i)%u_minus,faces(i)%u_plus,faces(i)%dn_minus,faces(i)%dn_plus, &
          faces(i)%h_normal,faces(i)%face_weight,penalty_factor,candidate(i),info)
        if(info /= 0)local_ok=.false.
      end if
    end do
    call collective_and(local_ok,communicator,ok)
    if (.not.ok) then
      message='nodal SIPG: canonical face evaluation failed collectively'
      deallocate(candidate)
      return
    end if
    call distribute_canonical_actions_mpi(faces,communicator,candidate,ok)
    if(.not.ok)then
      message='nodal SIPG: every physical face requires exactly one canonical owner'
      deallocate(candidate)
      return
    endif

    local_sums=0.0d0
    local_energy_imag=0.0d0
    local_balance=0.0d0
    do i=1,size(faces)
      if (.not.faces(i)%canonical_owner .or. .not.faces(i)%physical_face .or. &
          faces(i)%auxiliary_periodic_wrap) cycle
      energy=conjg(faces(i)%u_minus)*candidate(i)%total_value(1)+ &
        conjg(faces(i)%u_plus)*candidate(i)%total_value(2)+ &
        conjg(faces(i)%dn_minus)*candidate(i)%total_normal(1)+ &
        conjg(faces(i)%dn_plus)*candidate(i)%total_normal(2)
      local_energy_imag=max(local_energy_imag,abs(aimag(energy)))
      call evaluate_dg_nodal_sipg_face(faces(i)%u_plus,faces(i)%u_minus,-faces(i)%dn_plus,-faces(i)%dn_minus, &
        faces(i)%h_normal,faces(i)%face_weight,penalty_factor,reverse_action,info)
      local_balance=max(local_balance,maxval(abs(reverse_action%total_value-candidate(i)%total_value(2:1:-1))), &
        maxval(abs(reverse_action%total_normal+candidate(i)%total_normal(2:1:-1))))
      local_sums(1)=local_sums(1)+faces(i)%face_weight*abs(faces(i)%u_minus-faces(i)%u_plus)**2
      local_sums(2)=local_sums(2)+faces(i)%face_weight*penalty_factor/faces(i)%h_normal* &
        abs(faces(i)%u_minus-faces(i)%u_plus)**2
      local_sums(3)=local_sums(3)+1.0d0
    end do
#ifdef USE_MPI
    call MPI_Allreduce(local_energy_imag,candidate_diagnostics%hermiticity_defect,1, &
      MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok)then;ok=.false.;message='nodal SIPG: Hermiticity reduction failed collectively'; &
      deallocate(candidate);return;endif
    call MPI_Allreduce(local_balance,candidate_diagnostics%internal_cancellation_defect,1, &
      MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok)then;ok=.false.;message='nodal SIPG: balance reduction failed collectively'; &
      deallocate(candidate);return;endif
    call MPI_Allreduce(local_sums,global_sums,3,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok)then;ok=.false.;message='nodal SIPG: diagnostic reduction failed collectively'; &
      deallocate(candidate);return;endif
#else
    candidate_diagnostics%hermiticity_defect=local_energy_imag
    candidate_diagnostics%internal_cancellation_defect=local_balance
    global_sums=local_sums
    if(communicator < 0) then; ok=.false.; message='nodal SIPG: invalid serial communicator'; &
      deallocate(candidate); return; endif
#endif
    candidate_diagnostics%jump_norm=sqrt(global_sums(1))
    candidate_diagnostics%penalty_energy=global_sums(2)
    candidate_diagnostics%ownership_count=nint(global_sums(3))
    candidate_diagnostics%trace_epoch=expected_epoch
    if(allocated(action))deallocate(action)
    call move_alloc(candidate,action)
    diagnostics=candidate_diagnostics
    ok=.true.
    message=''
  end subroutine apply_dg_nodal_sipg_faces_mpi

  pure logical function valid_face(face,expected_epoch) result(ok)
    type(s_dg_nodal_sipg_face), intent(in) :: face
    integer(8), intent(in) :: expected_epoch
    real(8) :: normal_norm
    integer :: direction
    normal_norm=sqrt(sum(face%normal**2))
    ok=face%global_face_id > 0_8 .and. face%fragment_minus > 0 .and. face%fragment_plus > 0 .and. &
      face%axis >= 1 .and. face%axis <= 3 .and. face%h_normal > 0.0d0 .and. face%face_weight >= 0.0d0 .and. &
      ieee_is_finite(face%h_normal) .and. ieee_is_finite(face%face_weight) .and. &
      ieee_is_finite(normal_norm) .and. abs(normal_norm-1.0d0) <= 1.0d-12 .and. &
      finite_complex(face%u_minus) .and. finite_complex(face%u_plus) .and. &
      finite_complex(face%dn_minus) .and. finite_complex(face%dn_plus)
    if(ok)then
      do direction=1,3
        if(direction==face%axis)then
          ok=ok.and.abs(abs(face%normal(direction))-1.0d0)<=1.0d-12
        else
          ok=ok.and.abs(face%normal(direction))<=1.0d-12
        endif
      enddo
    endif
    if(face%physical_face)ok=ok .and. .not.face%auxiliary_periodic_wrap .and. face%trace_epoch==expected_epoch
    if(face%physical_supercell_periodic)then
      ok=ok.and.face%physical_face.and.abs(face%periodic_shift(face%axis))==1
      do direction=1,3
        if(direction/=face%axis)ok=ok.and.face%periodic_shift(direction)==0
      enddo
    else if(face%physical_face)then
      ok=ok.and.all(face%periodic_shift==0)
    endif
  end function valid_face

  subroutine distribute_canonical_actions_mpi(faces,communicator,actions,ok)
    type(s_dg_nodal_sipg_face),intent(in)::faces(:)
    integer, intent(in) :: communicator
    type(s_dg_nodal_sipg_action),intent(inout)::actions(:)
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer,parameter::action_width=10,metadata_integer_width=8,metadata_real_width=5,trace_width=4
    integer :: ierr,nproc,i,j,destination,nreceive,nentry,group_start,group_end,owner_position
    integer, allocatable :: send_count(:),receive_count(:),send_displacement(:),receive_displacement(:),cursor(:)
    integer, allocatable :: send_owner(:),receive_owner(:),send_index(:),receive_index(:),order(:)
    integer, allocatable :: complex_send_count(:),complex_receive_count(:),complex_send_displacement(:), &
      complex_receive_displacement(:),returned_index(:)
    integer(8), allocatable :: send_ids(:),receive_ids(:),send_metadata_integer(:),receive_metadata_integer(:)
    real(8),allocatable::send_metadata_real(:),receive_metadata_real(:)
    complex(8),allocatable::send_trace(:),receive_trace(:)
    complex(8),allocatable::send_action(:),receive_action(:),response_action(:),returned_action(:)
    logical :: local_ok,stage_ok
    call MPI_Comm_size(communicator,nproc,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok)then;ok=.false.;return;endif
    allocate(send_count(nproc),receive_count(nproc),send_displacement(nproc),receive_displacement(nproc))
    send_count=0
    do i=1,size(faces)
      if(.not.faces(i)%physical_face.or.faces(i)%auxiliary_periodic_wrap)cycle
      destination=int(modulo(faces(i)%global_face_id-1_8,int(nproc,8)))+1
      send_count(destination)=send_count(destination)+1
    end do
    call MPI_Alltoall(send_count,1,MPI_INTEGER,receive_count,1,MPI_INTEGER,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok)
    if(.not.stage_ok)then;ok=.false.;return;endif
    send_displacement(1)=0;receive_displacement(1)=0
    do i=2,nproc
      send_displacement(i)=send_displacement(i-1)+send_count(i-1)
      receive_displacement(i)=receive_displacement(i-1)+receive_count(i-1)
    end do
    nentry=sum(send_count)
    allocate(send_ids(nentry),send_owner(nentry),send_index(nentry),send_action(action_width*nentry), &
      send_metadata_integer(metadata_integer_width*nentry),send_metadata_real(metadata_real_width*nentry), &
      send_trace(trace_width*nentry),cursor(nproc))
    cursor=send_displacement+1
    do i=1,size(faces)
      if(.not.faces(i)%physical_face.or.faces(i)%auxiliary_periodic_wrap)cycle
      destination=int(modulo(faces(i)%global_face_id-1_8,int(nproc,8)))+1
      j=cursor(destination)
      send_ids(j)=faces(i)%global_face_id
      send_owner(j)=merge(1,0,faces(i)%canonical_owner)
      send_index(j)=i
      send_metadata_integer(metadata_integer_width*(j-1)+1:metadata_integer_width*j)= &
        [int(faces(i)%fragment_minus,8),int(faces(i)%fragment_plus,8),int(faces(i)%axis,8), &
         int(faces(i)%periodic_shift,8),merge(1_8,0_8,faces(i)%physical_supercell_periodic),faces(i)%trace_epoch]
      send_metadata_real(metadata_real_width*(j-1)+1:metadata_real_width*j)= &
        [faces(i)%normal,faces(i)%h_normal,faces(i)%face_weight]
      send_trace(trace_width*(j-1)+1:trace_width*j)= &
        [faces(i)%u_minus,faces(i)%u_plus,faces(i)%dn_minus,faces(i)%dn_plus]
      call pack_action(actions(i),send_action(action_width*(j-1)+1:action_width*j))
      cursor(destination)=j+1
    end do
    nreceive=sum(receive_count)
    allocate(receive_ids(nreceive),receive_owner(nreceive),receive_index(nreceive),order(nreceive))
    call MPI_Alltoallv(send_ids,send_count,send_displacement,MPI_INTEGER8,receive_ids,receive_count, &
      receive_displacement,MPI_INTEGER8,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    call MPI_Alltoallv(send_owner,send_count,send_displacement,MPI_INTEGER,receive_owner,receive_count, &
      receive_displacement,MPI_INTEGER,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    call MPI_Alltoallv(send_index,send_count,send_displacement,MPI_INTEGER,receive_index,receive_count, &
      receive_displacement,MPI_INTEGER,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    allocate(complex_send_count(nproc),complex_receive_count(nproc),complex_send_displacement(nproc), &
      complex_receive_displacement(nproc))
    allocate(receive_metadata_integer(metadata_integer_width*nreceive), &
      receive_metadata_real(metadata_real_width*nreceive),receive_trace(trace_width*nreceive))
    complex_send_count=metadata_integer_width*send_count
    complex_receive_count=metadata_integer_width*receive_count
    complex_send_displacement=metadata_integer_width*send_displacement
    complex_receive_displacement=metadata_integer_width*receive_displacement
    call MPI_Alltoallv(send_metadata_integer,complex_send_count,complex_send_displacement,MPI_INTEGER8, &
      receive_metadata_integer,complex_receive_count,complex_receive_displacement,MPI_INTEGER8,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    complex_send_count=metadata_real_width*send_count
    complex_receive_count=metadata_real_width*receive_count
    complex_send_displacement=metadata_real_width*send_displacement
    complex_receive_displacement=metadata_real_width*receive_displacement
    call MPI_Alltoallv(send_metadata_real,complex_send_count,complex_send_displacement,MPI_DOUBLE_PRECISION, &
      receive_metadata_real,complex_receive_count,complex_receive_displacement,MPI_DOUBLE_PRECISION,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    complex_send_count=trace_width*send_count
    complex_receive_count=trace_width*receive_count
    complex_send_displacement=trace_width*send_displacement
    complex_receive_displacement=trace_width*receive_displacement
    call MPI_Alltoallv(send_trace,complex_send_count,complex_send_displacement,MPI_DOUBLE_COMPLEX, &
      receive_trace,complex_receive_count,complex_receive_displacement,MPI_DOUBLE_COMPLEX,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    complex_send_count=action_width*send_count;complex_receive_count=action_width*receive_count
    complex_send_displacement=action_width*send_displacement
    complex_receive_displacement=action_width*receive_displacement
    allocate(receive_action(action_width*nreceive),response_action(action_width*nreceive))
    call MPI_Alltoallv(send_action,complex_send_count,complex_send_displacement,MPI_DOUBLE_COMPLEX, &
      receive_action,complex_receive_count,complex_receive_displacement,MPI_DOUBLE_COMPLEX,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    order=[(i,i=1,nreceive)]
    call heap_sort_face_order(receive_ids,order)
    response_action=(0.0d0,0.0d0);local_ok=.true.;group_start=1
    do while(group_start<=nreceive)
      group_end=group_start
      do while(group_end<nreceive)
        if(receive_ids(order(group_end+1))/=receive_ids(order(group_start)))exit
        group_end=group_end+1
      enddo
      owner_position=0
      do i=group_start,group_end
        if(receive_owner(order(i))==1)then
          if(owner_position/=0)local_ok=.false.
          owner_position=order(i)
        endif
      enddo
      if(owner_position==0)local_ok=.false.
      if(owner_position/=0)then
        do i=group_start,group_end
          j=order(i)
          if(.not.same_canonical_face_metadata(j,owner_position,receive_metadata_integer, &
            receive_metadata_real,receive_trace))local_ok=.false.
          response_action(action_width*(j-1)+1:action_width*j)= &
            receive_action(action_width*(owner_position-1)+1:action_width*owner_position)
        enddo
      endif
      group_start=group_end+1
    enddo
    call collective_and(local_ok,communicator,ok)
    if(.not.ok)return
    allocate(returned_index(nentry),returned_action(action_width*nentry))
    call MPI_Alltoallv(receive_index,receive_count,receive_displacement,MPI_INTEGER,returned_index,send_count, &
      send_displacement,MPI_INTEGER,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    call MPI_Alltoallv(response_action,complex_receive_count,complex_receive_displacement,MPI_DOUBLE_COMPLEX, &
      returned_action,complex_send_count,complex_send_displacement,MPI_DOUBLE_COMPLEX,communicator,ierr)
    call mpi_stage_agrees(ierr==MPI_SUCCESS,communicator,stage_ok);if(.not.stage_ok)then;ok=.false.;return;endif
    do i=1,nentry
      call unpack_action(returned_action(action_width*(i-1)+1:action_width*i),actions(returned_index(i)))
    enddo
#else
    integer::i,j,owners
    ok=communicator>=0
    do i=1,size(faces)
      if(.not.faces(i)%physical_face.or.faces(i)%auxiliary_periodic_wrap)cycle
      owners=0
      do j=1,size(faces)
        if(faces(j)%physical_face.and..not.faces(j)%auxiliary_periodic_wrap.and. &
          faces(j)%global_face_id==faces(i)%global_face_id.and.faces(j)%canonical_owner)owners=owners+1
      enddo
      if(owners/=1)ok=.false.
    enddo
#endif
  end subroutine distribute_canonical_actions_mpi

  subroutine collective_and(local_ok,communicator,global_ok)
    logical,intent(in)::local_ok
    integer,intent(in)::communicator
    logical,intent(out)::global_ok
#ifdef USE_MPI
    integer::ierr
    call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS)global_ok=.false.
#else
    global_ok=local_ok.and.communicator>=0
#endif
  end subroutine collective_and

  pure subroutine heap_sort_face_order(ids,order)
    integer(8),intent(in)::ids(:)
    integer,intent(inout)::order(:)
    integer::root,tail,child,saved
    if(size(order)<2)return
    do root=size(order)/2,1,-1
      saved=order(root);tail=root
      do while(2*tail<=size(order))
        child=2*tail
        if(child<size(order))then
          if(ids(order(child))<ids(order(child+1)))child=child+1
        endif
        if(ids(saved)>=ids(order(child)))exit
        order(tail)=order(child);tail=child
      enddo
      order(tail)=saved
    enddo
    do root=size(order),2,-1
      saved=order(root);order(root)=order(1);tail=1
      do while(2*tail<root)
        child=2*tail
        if(child<root-1)then
          if(ids(order(child))<ids(order(child+1)))child=child+1
        endif
        if(ids(saved)>=ids(order(child)))exit
        order(tail)=order(child);tail=child
      enddo
      order(tail)=saved
    enddo
  end subroutine heap_sort_face_order

  pure subroutine pack_action(action,packed)
    type(s_dg_nodal_sipg_action),intent(in)::action
    complex(8),intent(out)::packed(10)
    packed=[action%consistency_value,action%symmetry_normal,action%penalty_value, &
      action%total_value,action%total_normal]
  end subroutine pack_action

  pure subroutine unpack_action(packed,action)
    complex(8),intent(in)::packed(10)
    type(s_dg_nodal_sipg_action),intent(out)::action
    action%consistency_value=packed(1:2)
    action%symmetry_normal=packed(3:4)
    action%penalty_value=packed(5:6)
    action%total_value=packed(7:8)
    action%total_normal=packed(9:10)
  end subroutine unpack_action

  pure logical function same_canonical_face_metadata(left,right,integer_metadata,real_metadata,traces) result(same)
    integer,intent(in)::left,right
    integer(8),intent(in)::integer_metadata(:)
    real(8),intent(in)::real_metadata(:)
    complex(8),intent(in)::traces(:)
    real(8)::real_scale,trace_scale
    same=all(integer_metadata(8*(left-1)+1:8*left)==integer_metadata(8*(right-1)+1:8*right))
    real_scale=max(1.0d0,maxval(abs(real_metadata(5*(right-1)+1:5*right))))
    same=same.and.maxval(abs(real_metadata(5*(left-1)+1:5*left)- &
      real_metadata(5*(right-1)+1:5*right)))<=1.0d-13*real_scale
    trace_scale=max(1.0d0,maxval(abs(traces(4*(right-1)+1:4*right))))
    same=same.and.maxval(abs(traces(4*(left-1)+1:4*left)- &
      traces(4*(right-1)+1:4*right)))<=1.0d-13*trace_scale
  end function same_canonical_face_metadata

#ifdef USE_MPI
  subroutine mpi_stage_agrees(local_ok,communicator,global_ok)
    logical,intent(in)::local_ok
    integer,intent(in)::communicator
    logical,intent(out)::global_ok
    integer::ierr
    call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS)global_ok=.false.
  end subroutine mpi_stage_agrees
#endif

  pure logical function finite_complex(value) result(ok)
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function finite_complex
end module dg_nodal_sipg
