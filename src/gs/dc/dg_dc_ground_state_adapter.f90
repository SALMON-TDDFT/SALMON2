!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_dc_ground_state_adapter
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dg_nodal_state, only: s_dg_nodal_common_state,dg_nodal_face_slot,allocate_dg_nodal_face_buffers
  use dg_nodal_sipg, only: s_dg_nodal_sipg_action,evaluate_dg_nodal_sipg_face
  use dg_dc_ground_state, only: s_dg_dc_gs_diagnostics
  use dg_dc_buffer_core_faces, only: s_dg_dc_buffer_core_face,build_dg_dc_buffer_core_faces, &
    exchange_dg_dc_buffer_core_state_window,validate_dg_dc_buffer_core_faces, &
    assemble_dg_dc_local_buffer_state_window,project_dg_dc_buffer_core_face
  implicit none
  private
  real(8), parameter :: dg_dc_maximum_candidate_workspace_bytes=8d0*1024d0**3
  public :: compose_dg_dc_hamiltonian
  public :: reconstruct_dg_dc_core_density
  public :: mix_dg_dc_density_transaction
  public :: rollback_dg_dc_density_potential
  public :: lift_dg_dc_sipg_face
  public :: expand_dg_dc_global_candidate_axis
  public :: build_dg_dc_one_sided_derivative_weights
  public :: evaluate_dg_dc_local_sipg_face
  public :: exchange_dg_dc_face_traces_mpi
  public :: find_dg_dc_periodic_face_neighbor
  public :: initialize_dg_dc_physical_faces,apply_dg_dc_sipg_operator_mpi
  public :: build_dg_dc_interior_volume_action
  public :: measure_dg_dc_face_action_balance
  public :: measure_dg_dc_cross_hermiticity_mpi
  public :: validate_dg_dc_candidate_memory
  public :: execute_dg_dc_production_iteration
  public :: prepare_dg_dc_buffer_core_projectors

  abstract interface
    subroutine dg_dc_restore_callback(state,density,communicator,ok,message)
      import :: s_dg_nodal_common_state
      type(s_dg_nodal_common_state), intent(inout) :: state
      real(8), intent(inout) :: density(:,:,:,:)
      integer, intent(in) :: communicator
      logical, intent(out) :: ok
      character(*), intent(out) :: message
    end subroutine dg_dc_restore_callback

    subroutine dg_dc_solve_callback(state,lambda,communicator,diagnostics,ok,message)
      import :: s_dg_nodal_common_state,s_dg_dc_gs_diagnostics
      type(s_dg_nodal_common_state), intent(inout) :: state
      real(8), intent(in) :: lambda
      integer, intent(in) :: communicator
      type(s_dg_dc_gs_diagnostics), intent(inout) :: diagnostics
      logical, intent(out) :: ok
      character(*), intent(out) :: message
    end subroutine dg_dc_solve_callback

    subroutine dg_dc_update_callback(state,density,density_mix,unmixed,communicator,diagnostics,ok,message)
      import :: s_dg_nodal_common_state,s_dg_dc_gs_diagnostics
      type(s_dg_nodal_common_state), intent(inout) :: state
      real(8), intent(inout) :: density(:,:,:,:)
      real(8), intent(in) :: density_mix
      logical, intent(in) :: unmixed
      integer, intent(in) :: communicator
      type(s_dg_dc_gs_diagnostics), intent(inout) :: diagnostics
      logical, intent(out) :: ok
      character(*), intent(out) :: message
    end subroutine dg_dc_update_callback
  end interface

contains

  subroutine prepare_dg_dc_buffer_core_projectors(local_fragment,fragment_origins,fragment_sizes,&
      global_size,local_core_origin,local_core_size,buffer_width,generation,operator_fingerprint,&
      core_physical_grid_ids,core_values,core_gradients,owned_state_ids,configured_states,&
      buffer_face_values,buffer_face_weights,rank_tolerance,minimum_retained_rank,&
      maximum_projection_residual,maximum_escape_norm,maximum_projection_mismatch,&
      communicator,faces,ok,message)
    integer,intent(in) :: local_fragment,fragment_origins(:,:),fragment_sizes(:,:),global_size(3)
    integer,intent(in) :: local_core_origin(3),local_core_size(3),buffer_width,generation
    integer(8),intent(in) :: operator_fingerprint,core_physical_grid_ids(:,:,:)
    real(8),intent(in) :: core_values(:,:,:,:),core_gradients(:,:,:,:,:)
    integer,intent(in) :: owned_state_ids(:),configured_states
    real(8),intent(in) :: buffer_face_values(:,:,:),buffer_face_weights(:,:)
    real(8),intent(in) :: rank_tolerance
    integer,intent(in) :: minimum_retained_rank
    real(8),intent(in) :: maximum_projection_residual,maximum_escape_norm,maximum_projection_mismatch
    integer,intent(in) :: communicator
    type(s_dg_dc_buffer_core_face),allocatable,intent(out) :: faces(:)
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: iface
    logical :: local_ok

    call build_dg_dc_buffer_core_faces(local_fragment,fragment_origins,fragment_sizes,global_size,&
      buffer_width,generation,communicator,faces,local_ok,message,local_core_origin,local_core_size)
    call dg_dc_collective_and(local_ok,communicator,ok)
    if(.not.ok)return
    call validate_dg_dc_buffer_core_faces(faces,fragment_origins,fragment_sizes,global_size,&
      communicator,local_ok,message)
    call dg_dc_collective_and(local_ok,communicator,ok)
    if(.not.ok)return
    call exchange_dg_dc_buffer_core_state_window(faces,fragment_origins,fragment_sizes,global_size,&
      local_core_origin,local_core_size,&
      core_physical_grid_ids,core_values,core_gradients,owned_state_ids,configured_states,&
      communicator,local_ok,message)
    call dg_dc_collective_and(local_ok,communicator,ok)
    if(.not.ok)return
    call assemble_dg_dc_local_buffer_state_window(faces,buffer_face_values,owned_state_ids,&
      configured_states,communicator,local_ok,message)
    call dg_dc_collective_and(local_ok,communicator,ok)
    if(.not.ok)return
    local_ok=size(buffer_face_values,2)==size(owned_state_ids).and.size(buffer_face_values,3)==6.and.&
      size(buffer_face_weights,2)==6.and.size(buffer_face_values,1)==size(buffer_face_weights,1)
    do iface=1,6
      if(.not.local_ok)exit
      if(faces(iface)%point_count==0)cycle
      if(faces(iface)%point_count>size(buffer_face_values,1))then
        local_ok=.false.;exit
      endif
      call project_dg_dc_buffer_core_face(faces(iface),&
        faces(iface)%local_buffer_values,&
        buffer_face_weights(1:faces(iface)%point_count,iface),rank_tolerance,minimum_retained_rank,&
        maximum_projection_residual,maximum_escape_norm,maximum_projection_mismatch,generation,&
        operator_fingerprint,local_ok,message)
      if(.not.local_ok)exit
    enddo
    call dg_dc_collective_and(local_ok,communicator,ok)
    if(ok)then
      message=''
    else if(len_trim(message)==0)then
      message='DG DC GS adapter: collective buffer/core projector preparation failed'
    endif
  end subroutine prepare_dg_dc_buffer_core_projectors

  subroutine execute_dg_dc_production_iteration(state,density,lambda,density_mix,reset_mixing,unmixed, &
      restore_potential,solve_complete_h,update_density_potential,communicator,diagnostics,ok,message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    real(8), intent(inout) :: density(:,:,:,:)
    real(8), intent(in) :: lambda,density_mix
    logical, intent(in) :: reset_mixing,unmixed
    procedure(dg_dc_restore_callback) :: restore_potential
    procedure(dg_dc_solve_callback) :: solve_complete_h
    procedure(dg_dc_update_callback) :: update_density_potential
    integer, intent(in) :: communicator
    type(s_dg_dc_gs_diagnostics), intent(out) :: diagnostics
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    logical :: stage_ok
    diagnostics=s_dg_dc_gs_diagnostics()
    if(reset_mixing) then
      call restore_potential(state,density,communicator,stage_ok,message)
      call dg_dc_collective_and(stage_ok,communicator,ok)
      if(.not.ok) return
    end if
    call solve_complete_h(state,lambda,communicator,diagnostics,stage_ok,message)
    call dg_dc_collective_and(stage_ok,communicator,ok)
    if(.not.ok) return
    call update_density_potential(state,density,density_mix,unmixed,communicator,diagnostics,stage_ok,message)
    call dg_dc_collective_and(stage_ok,communicator,ok)
    if(ok) message=''
  end subroutine execute_dg_dc_production_iteration

  subroutine compose_dg_dc_hamiltonian(hpsi_volume_nonlocal,hpsi_sipg,lambda,hpsi,ok,message)
    complex(8), intent(in) :: hpsi_volume_nonlocal(:,:,:,:,:),hpsi_sipg(:,:,:,:,:)
    real(8), intent(in) :: lambda
    complex(8), intent(out) :: hpsi(:,:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    ok=lambda>=0d0 .and. lambda<=1d0 .and. all(shape(hpsi)==shape(hpsi_volume_nonlocal)) .and. &
      all(shape(hpsi)==shape(hpsi_sipg)) .and. all(ieee_is_finite(real(hpsi_volume_nonlocal,8))) .and. &
      all(ieee_is_finite(aimag(hpsi_volume_nonlocal))) .and. all(ieee_is_finite(real(hpsi_sipg,8))) .and. &
      all(ieee_is_finite(aimag(hpsi_sipg)))
    if(.not.ok) then
      hpsi=(0d0,0d0)
      message='DG DC GS adapter: invalid Hamiltonian action'
      return
    end if
    hpsi=hpsi_volume_nonlocal+lambda*hpsi_sipg
    ok=all(ieee_is_finite(real(hpsi,8))) .and. all(ieee_is_finite(aimag(hpsi)))
    if(ok) then
      message=''
    else
      hpsi=(0d0,0d0)
      message='DG DC GS adapter: nonfinite composed Hamiltonian'
    end if
  end subroutine compose_dg_dc_hamiltonian

  subroutine reconstruct_dg_dc_core_density(state,valid_candidate_count,density,ok,message)
    type(s_dg_nodal_common_state), intent(in) :: state
    integer, intent(in) :: valid_candidate_count
    real(8), intent(out) :: density(:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: istate,ispin
    ok=state%initialized .and. allocated(state%psi_core) .and. allocated(state%occupations) .and. &
      valid_candidate_count>=1 .and. valid_candidate_count<=state%nstate .and. &
      all(shape(density)==[state%core_size,state%nspin])
    if(.not.ok) then
      density=0d0
      message='DG DC GS adapter: invalid core-density layout'
      return
    end if
    density=0d0
    do ispin=1,state%nspin
      do istate=1,valid_candidate_count
        density(:,:,:,ispin)=density(:,:,:,ispin)+state%occupations(istate,ispin)* &
          abs(state%psi_core(:,:,:,istate,ispin))**2
      end do
    end do
    ok=all(ieee_is_finite(density)) .and. all(density>=0d0)
    if(ok) then
      message=''
    else
      density=0d0
      message='DG DC GS adapter: nonfinite core density'
    end if
  end subroutine reconstruct_dg_dc_core_density

  subroutine mix_dg_dc_density_transaction(accepted_density,raw_density,mix_rate,unmixed,mixed_density,ok,message)
    real(8), intent(in) :: accepted_density(:,:,:,:),raw_density(:,:,:,:)
    real(8), intent(in) :: mix_rate
    logical, intent(in) :: unmixed
    real(8), intent(out) :: mixed_density(:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: effective_mix
    ok=all(shape(accepted_density)==shape(raw_density)) .and. &
      all(shape(mixed_density)==shape(raw_density)) .and. mix_rate>0d0 .and. mix_rate<=1d0 .and. &
      all(ieee_is_finite(accepted_density)) .and. all(ieee_is_finite(raw_density))
    if(.not.ok) then
      mixed_density=0d0
      message='DG DC GS adapter: invalid density mixing transaction'
      return
    end if
    effective_mix=merge(1d0,mix_rate,unmixed)
    mixed_density=(1d0-effective_mix)*accepted_density+effective_mix*raw_density
    ok=all(ieee_is_finite(mixed_density))
    if(ok) then
      message=''
    else
      mixed_density=0d0
      message='DG DC GS adapter: nonfinite mixed density'
    end if
  end subroutine mix_dg_dc_density_transaction

  subroutine rollback_dg_dc_density_potential(accepted_density,accepted_potential,density,potential)
    real(8), intent(in) :: accepted_density(:,:,:,:),accepted_potential(:,:,:,:)
    real(8), intent(out) :: density(:,:,:,:),potential(:,:,:,:)
    density=accepted_density
    potential=accepted_potential
  end subroutine rollback_dg_dc_density_potential

  subroutine lift_dg_dc_sipg_face(axis,side,derivative_weights,face_value,face_normal,hpsi,ok,message)
    integer, intent(in) :: axis,side
    real(8), intent(in) :: derivative_weights(:)
    complex(8), intent(in) :: face_value(:,:,:,:),face_normal(:,:,:,:)
    complex(8), intent(inout) :: hpsi(:,:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: i,j,k,istate,ispin,index_normal,nnormal
    integer :: expected_shape(4)
    nnormal=size(derivative_weights)
    select case(axis)
    case(1)
      expected_shape=[size(hpsi,2),size(hpsi,3),size(hpsi,4),size(hpsi,5)]
      ok=nnormal<=size(hpsi,1)
    case(2)
      expected_shape=[size(hpsi,1),size(hpsi,3),size(hpsi,4),size(hpsi,5)]
      ok=nnormal<=size(hpsi,2)
    case(3)
      expected_shape=[size(hpsi,1),size(hpsi,2),size(hpsi,4),size(hpsi,5)]
      ok=nnormal<=size(hpsi,3)
    case default
      expected_shape=0
      ok=.false.
    end select
    ok=ok .and. abs(side)==1 .and. nnormal>=1 .and. all(shape(face_value)==expected_shape) .and. &
      all(shape(face_normal)==expected_shape) .and. all(ieee_is_finite(derivative_weights)) .and. &
      all(ieee_is_finite(real(face_value,8))) .and. all(ieee_is_finite(aimag(face_value))) .and. &
      all(ieee_is_finite(real(face_normal,8))) .and. all(ieee_is_finite(aimag(face_normal)))
    if(.not.ok) then
      message='DG DC GS adapter: invalid SIPG face lift'
      return
    end if
    do ispin=1,size(hpsi,5)
      do istate=1,size(hpsi,4)
        select case(axis)
        case(1)
          do k=1,size(hpsi,3)
            do j=1,size(hpsi,2)
              i=merge(1,size(hpsi,1),side<0)
              hpsi(i,j,k,istate,ispin)=hpsi(i,j,k,istate,ispin)+face_value(j,k,istate,ispin)
              do index_normal=1,nnormal
                i=merge(index_normal,size(hpsi,1)-index_normal+1,side<0)
                hpsi(i,j,k,istate,ispin)=hpsi(i,j,k,istate,ispin)+ &
                  derivative_weights(index_normal)*face_normal(j,k,istate,ispin)
              end do
            end do
          end do
        case(2)
          do k=1,size(hpsi,3)
            do i=1,size(hpsi,1)
              j=merge(1,size(hpsi,2),side<0)
              hpsi(i,j,k,istate,ispin)=hpsi(i,j,k,istate,ispin)+face_value(i,k,istate,ispin)
              do index_normal=1,nnormal
                j=merge(index_normal,size(hpsi,2)-index_normal+1,side<0)
                hpsi(i,j,k,istate,ispin)=hpsi(i,j,k,istate,ispin)+ &
                  derivative_weights(index_normal)*face_normal(i,k,istate,ispin)
              end do
            end do
          end do
        case(3)
          do j=1,size(hpsi,2)
            do i=1,size(hpsi,1)
              k=merge(1,size(hpsi,3),side<0)
              hpsi(i,j,k,istate,ispin)=hpsi(i,j,k,istate,ispin)+face_value(i,j,istate,ispin)
              do index_normal=1,nnormal
                k=merge(index_normal,size(hpsi,3)-index_normal+1,side<0)
                hpsi(i,j,k,istate,ispin)=hpsi(i,j,k,istate,ispin)+ &
                  derivative_weights(index_normal)*face_normal(i,j,istate,ispin)
              end do
            end do
          end do
        end select
      end do
    end do
    ok=all(ieee_is_finite(real(hpsi,8))) .and. all(ieee_is_finite(aimag(hpsi)))
    if(ok) then
      message=''
    else
      message='DG DC GS adapter: nonfinite lifted SIPG action'
    end if
  end subroutine lift_dg_dc_sipg_face

  subroutine expand_dg_dc_global_candidate_axis(state,local_candidate_count,communicator,ok,message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    integer, intent(in) :: local_candidate_count,communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: psi_global(:,:,:,:,:),hpsi_global(:,:,:,:,:)
    real(8), allocatable :: occupations_global(:,:)
    integer :: global_candidate_count,candidate_offset,allocation_status
    real(8) :: estimated_bytes
    logical :: memory_ok
#ifdef USE_MPI
    integer :: ierr,ierr_exscan,ierr_rank,rank,local_valid,global_valid
    local_valid=merge(1,0,state%initialized .and. local_candidate_count>=1 .and. &
      local_candidate_count<=state%nstate)
    call MPI_Allreduce(local_valid,global_valid,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    ok=ierr==MPI_SUCCESS .and. global_valid==1
    if(.not.ok) then
      message='DG DC GS adapter: invalid local candidate count collectively'
      return
    end if
    call MPI_Allreduce(local_candidate_count,global_candidate_count,1,MPI_INTEGER,MPI_SUM,communicator,ierr)
    ok=ierr==MPI_SUCCESS .and. global_candidate_count>=local_candidate_count
    if(.not.ok) then
      message='DG DC GS adapter: global candidate count reduction failed'
      return
    end if
    candidate_offset=0
    rank=-1
    call MPI_Exscan(local_candidate_count,candidate_offset,1,MPI_INTEGER,MPI_SUM,communicator,ierr_exscan)
    call MPI_Comm_rank(communicator,rank,ierr_rank)
    ok=ierr_exscan==MPI_SUCCESS .and. ierr_rank==MPI_SUCCESS
    if(ok .and. rank==0) candidate_offset=0
    call dg_dc_collective_and(ok,communicator,memory_ok)
    ok=memory_ok
#else
    ok=state%initialized .and. local_candidate_count>=1 .and. local_candidate_count<=state%nstate .and. &
      communicator>=0
    global_candidate_count=local_candidate_count
    candidate_offset=0
#endif
    if(.not.ok) then
      message='DG DC GS adapter: candidate offset construction failed'
      return
    end if
    call validate_dg_dc_candidate_memory(state%core_size,global_candidate_count,state%nspin, &
      dg_dc_maximum_candidate_workspace_bytes,estimated_bytes,memory_ok)
    call dg_dc_collective_and(memory_ok,communicator,ok)
    if(.not.ok) then
      message='DG DC GS adapter: candidate workspace exceeds bounded per-rank memory'
      return
    end if
    allocate(psi_global(state%core_size(1),state%core_size(2),state%core_size(3), &
      global_candidate_count,state%nspin),stat=allocation_status)
    memory_ok=allocation_status==0
    if(memory_ok) then
      allocate(hpsi_global,mold=psi_global,stat=allocation_status)
      memory_ok=allocation_status==0
    end if
    if(memory_ok) then
      allocate(occupations_global(global_candidate_count,state%nspin),stat=allocation_status)
      memory_ok=allocation_status==0
    end if
    call dg_dc_collective_and(memory_ok,communicator,ok)
    if(.not.ok) then
      if(allocated(psi_global)) deallocate(psi_global)
      if(allocated(hpsi_global)) deallocate(hpsi_global)
      if(allocated(occupations_global)) deallocate(occupations_global)
      message='DG DC GS adapter: candidate workspace allocation failed collectively'
      return
    end if
    psi_global=(0d0,0d0)
    hpsi_global=(0d0,0d0)
    occupations_global=0d0
    psi_global(:,:,:,candidate_offset+1:candidate_offset+local_candidate_count,:)= &
      state%psi_core(:,:,:,1:local_candidate_count,:)
    hpsi_global(:,:,:,candidate_offset+1:candidate_offset+local_candidate_count,:)= &
      state%hpsi_core(:,:,:,1:local_candidate_count,:)
    occupations_global(candidate_offset+1:candidate_offset+local_candidate_count,:)= &
      state%occupations(1:local_candidate_count,:)
    call move_alloc(psi_global,state%psi_core)
    call move_alloc(hpsi_global,state%hpsi_core)
    call move_alloc(occupations_global,state%occupations)
    state%nstate=global_candidate_count
    message=''
  end subroutine expand_dg_dc_global_candidate_axis

  subroutine build_dg_dc_one_sided_derivative_weights(grid_spacing,side,weights,ok,message)
    real(8), intent(in) :: grid_spacing
    integer, intent(in) :: side
    real(8), intent(out) :: weights(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), allocatable :: matrix(:,:),rhs(:)
    real(8) :: pivot,factor
    integer :: n,row,column,pivot_row
    n=size(weights)
    ok=n>=2 .and. grid_spacing>0d0 .and. ieee_is_finite(grid_spacing) .and. abs(side)==1
    if(.not.ok) then
      weights=0d0
      message='DG DC GS adapter: invalid one-sided derivative stencil'
      return
    end if
    allocate(matrix(n,n),rhs(n))
    do row=1,n
      do column=1,n
        matrix(row,column)=(grid_spacing*dble(column-1))**(row-1)
      end do
    end do
    rhs=0d0
    rhs(2)=-1d0
    do column=1,n
      pivot_row=column-1+maxloc(abs(matrix(column:n,column)),dim=1)
      if(abs(matrix(pivot_row,column))<=tiny(1d0)) then
        ok=.false.; weights=0d0; message='DG DC GS adapter: singular derivative stencil'; return
      end if
      if(pivot_row/=column) then
        do row=column,n
          pivot=matrix(column,row)
          matrix(column,row)=matrix(pivot_row,row)
          matrix(pivot_row,row)=pivot
        end do
        pivot=rhs(column); rhs(column)=rhs(pivot_row); rhs(pivot_row)=pivot
      end if
      pivot=matrix(column,column)
      matrix(column,column:n)=matrix(column,column:n)/pivot
      rhs(column)=rhs(column)/pivot
      do row=1,n
        if(row==column) cycle
        factor=matrix(row,column)
        matrix(row,column:n)=matrix(row,column:n)-factor*matrix(column,column:n)
        rhs(row)=rhs(row)-factor*rhs(column)
      end do
    end do
    weights=rhs
    ok=all(ieee_is_finite(weights))
    if(ok) then
      message=''
    else
      weights=0d0
      message='DG DC GS adapter: nonfinite derivative weights'
    end if
  end subroutine build_dg_dc_one_sided_derivative_weights

  subroutine evaluate_dg_dc_local_sipg_face(local_value,local_outward_normal,neighbor_value, &
      neighbor_outward_normal,h_normal,face_weight,penalty_factor,side,local_face_value, &
      local_face_normal,ok,message)
    complex(8), intent(in) :: local_value(:,:,:,:),local_outward_normal(:,:,:,:)
    complex(8), intent(in) :: neighbor_value(:,:,:,:),neighbor_outward_normal(:,:,:,:)
    real(8), intent(in) :: h_normal,face_weight,penalty_factor
    integer, intent(in) :: side
    complex(8), intent(out) :: local_face_value(:,:,:,:),local_face_normal(:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_sipg_action) :: action
    integer :: i,j,istate,ispin,info
    ok=abs(side)==1 .and. all(shape(local_value)==shape(local_outward_normal)) .and. &
      all(shape(local_value)==shape(neighbor_value)) .and. &
      all(shape(local_value)==shape(neighbor_outward_normal)) .and. &
      all(shape(local_value)==shape(local_face_value)) .and. &
      all(shape(local_value)==shape(local_face_normal))
    if(.not.ok) then
      local_face_value=(0d0,0d0); local_face_normal=(0d0,0d0)
      message='DG DC GS adapter: invalid paired SIPG trace layout'
      return
    end if
    do ispin=1,size(local_value,4)
      do istate=1,size(local_value,3)
        do j=1,size(local_value,2)
          do i=1,size(local_value,1)
            if(side>0) then
              call evaluate_dg_nodal_sipg_face(local_value(i,j,istate,ispin), &
                neighbor_value(i,j,istate,ispin),local_outward_normal(i,j,istate,ispin), &
                -neighbor_outward_normal(i,j,istate,ispin),h_normal,face_weight,penalty_factor,action,info)
              local_face_value(i,j,istate,ispin)=action%total_value(1)
              local_face_normal(i,j,istate,ispin)=action%total_normal(1)
            else
              call evaluate_dg_nodal_sipg_face(neighbor_value(i,j,istate,ispin), &
                local_value(i,j,istate,ispin),neighbor_outward_normal(i,j,istate,ispin), &
                -local_outward_normal(i,j,istate,ispin),h_normal,face_weight,penalty_factor,action,info)
              local_face_value(i,j,istate,ispin)=action%total_value(2)
              local_face_normal(i,j,istate,ispin)=-action%total_normal(2)
            end if
            if(info/=0) ok=.false.
          end do
        end do
      end do
    end do
    if(ok) then
      message=''
    else
      local_face_value=(0d0,0d0); local_face_normal=(0d0,0d0)
      message='DG DC GS adapter: paired SIPG trace evaluation failed'
    end if
  end subroutine evaluate_dg_dc_local_sipg_face

  subroutine exchange_dg_dc_face_traces_mpi(local_value,local_normal,neighbor_rank,tag,communicator, &
      neighbor_value,neighbor_normal,ok,message)
    complex(8), intent(in) :: local_value(:,:,:,:),local_normal(:,:,:,:)
    integer, intent(in) :: neighbor_rank,tag,communicator
    complex(8), intent(out) :: neighbor_value(:,:,:,:),neighbor_normal(:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
#ifdef USE_MPI
    integer :: ierr_value,ierr_normal,ierr_agree,local_ok,global_ok
    ok=neighbor_rank>=0 .and. tag>=0 .and. all(shape(local_value)==shape(local_normal)) .and. &
      all(shape(local_value)==shape(neighbor_value)) .and. all(shape(local_value)==shape(neighbor_normal))
    local_ok=merge(1,0,ok)
    call MPI_Allreduce(local_ok,global_ok,1,MPI_INTEGER,MPI_MIN,communicator,ierr_agree)
    if(ierr_agree/=MPI_SUCCESS .or. global_ok/=1) then
      ok=.false.; neighbor_value=(0d0,0d0); neighbor_normal=(0d0,0d0)
      message='DG DC GS adapter: invalid face exchange layout collectively'
      return
    end if
    call MPI_Sendrecv(local_value,size(local_value),MPI_DOUBLE_COMPLEX,neighbor_rank,tag, &
      neighbor_value,size(neighbor_value),MPI_DOUBLE_COMPLEX,neighbor_rank,tag,communicator, &
      MPI_STATUS_IGNORE,ierr_value)
    call MPI_Sendrecv(local_normal,size(local_normal),MPI_DOUBLE_COMPLEX,neighbor_rank,tag+1, &
      neighbor_normal,size(neighbor_normal),MPI_DOUBLE_COMPLEX,neighbor_rank,tag+1,communicator, &
      MPI_STATUS_IGNORE,ierr_normal)
    local_ok=merge(1,0,ierr_value==MPI_SUCCESS .and. ierr_normal==MPI_SUCCESS)
    call MPI_Allreduce(local_ok,global_ok,1,MPI_INTEGER,MPI_MIN,communicator,ierr_agree)
    ok=ierr_agree==MPI_SUCCESS .and. global_ok==1
#else
    ok=neighbor_rank==0 .and. tag>=0 .and. communicator>=0 .and. &
      all(shape(local_value)==shape(local_normal)) .and. &
      all(shape(local_value)==shape(neighbor_value)) .and. all(shape(local_value)==shape(neighbor_normal))
    if(ok) then
      neighbor_value=local_value
      neighbor_normal=local_normal
    end if
#endif
    if(ok) then
      message=''
    else
      neighbor_value=(0d0,0d0); neighbor_normal=(0d0,0d0)
      message='DG DC GS adapter: face trace exchange failed collectively'
    end if
  end subroutine exchange_dg_dc_face_traces_mpi

  subroutine find_dg_dc_periodic_face_neighbor(fragment,axis,side,fragment_origins,fragment_sizes, &
      global_size,neighbor_fragment,ok,message)
    integer, intent(in) :: fragment,axis,side
    integer, intent(in) :: fragment_origins(:,:),fragment_sizes(:,:),global_size(3)
    integer, intent(out) :: neighbor_fragment
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: candidate,tangent,target
    logical :: matches
    ok=axis>=1 .and. axis<=3 .and. abs(side)==1 .and. fragment>=1 .and. &
      fragment<=size(fragment_origins,2) .and. all(shape(fragment_origins)==shape(fragment_sizes)) .and. &
      size(fragment_origins,1)==3 .and. all(global_size>0) .and. all(fragment_sizes>0)
    neighbor_fragment=0
    if(.not.ok) then
      message='DG DC GS adapter: invalid fragment face topology'
      return
    end if
    if(side>0) then
      target=modulo(fragment_origins(axis,fragment)+fragment_sizes(axis,fragment),global_size(axis))
    else
      target=fragment_origins(axis,fragment)
    end if
    do candidate=1,size(fragment_origins,2)
      matches=.true.
      do tangent=1,3
        if(tangent==axis) cycle
        matches=matches .and. fragment_origins(tangent,candidate)==fragment_origins(tangent,fragment) .and. &
          fragment_sizes(tangent,candidate)==fragment_sizes(tangent,fragment)
      end do
      if(side>0) then
        matches=matches .and. fragment_origins(axis,candidate)==target
      else
        matches=matches .and. modulo(fragment_origins(axis,candidate)+fragment_sizes(axis,candidate), &
          global_size(axis))==target
      end if
      if(matches) then
        if(neighbor_fragment/=0) then
          ok=.false.; message='DG DC GS adapter: ambiguous physical face neighbor'; return
        end if
        neighbor_fragment=candidate
      end if
    end do
    ok=neighbor_fragment/=0
    if(ok) then
      message=''
    else
      message='DG DC GS adapter: missing physical face neighbor'
    end if
  end subroutine find_dg_dc_periodic_face_neighbor

  subroutine initialize_dg_dc_physical_faces(state,fragment_origins,fragment_sizes,global_size, &
      communicator,ok,message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    integer, intent(in) :: fragment_origins(:,:),fragment_sizes(:,:),global_size(3),communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: axis,side,iface,neighbor
    logical :: face_ok,local_ok
#ifdef USE_MPI
    integer :: rank,nproc,ierr_rank,ierr_size
    rank=-1
    nproc=-1
    call MPI_Comm_rank(communicator,rank,ierr_rank)
    call MPI_Comm_size(communicator,nproc,ierr_size)
    local_ok=ierr_rank==MPI_SUCCESS .and. ierr_size==MPI_SUCCESS .and. &
      nproc==size(fragment_origins,2) .and. state%fragment>=1 .and. &
      state%fragment<=size(fragment_origins,2)
    if(local_ok) local_ok=rank==state%fragment-1 .and. &
      all(fragment_sizes(:,state%fragment)==state%core_size)
    call dg_dc_collective_and(local_ok,communicator,ok)
#else
    ok=communicator>=0 .and. size(fragment_origins,2)==1 .and. state%fragment==1
    if(ok) ok=all(fragment_sizes(:,state%fragment)==state%core_size)
#endif
    if(.not.ok) then
      message='DG DC GS adapter: fragment/rank topology is not one-to-one'
      return
    end if
    do axis=1,3
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        if(global_size(axis)==1) then
          neighbor=state%fragment
          face_ok=.true.
        else
          call find_dg_dc_periodic_face_neighbor(state%fragment,axis,side,fragment_origins, &
            fragment_sizes,global_size,neighbor,face_ok,message)
        end if
        call dg_dc_collective_and(face_ok,communicator,ok)
        if(.not.ok) then
          message='DG DC GS adapter: physical face neighbor construction failed collectively'
          return
        end if
        state%faces(iface)%neighbor_fragment=neighbor
        state%faces(iface)%neighbor_rank=neighbor-1
        call allocate_dg_nodal_face_buffers(state%faces(iface),state%core_size,1,state%nstate,state%nspin,face_ok)
        call dg_dc_collective_and(face_ok,communicator,ok)
        if(.not.ok) then
          message='DG DC GS adapter: physical face buffer allocation failed collectively'
          return
        end if
      end do
    end do
    ok=.true.; message=''
  end subroutine initialize_dg_dc_physical_faces

  subroutine apply_dg_dc_sipg_operator_mpi(state,fragment_origins,fragment_sizes,global_size,grid_spacing, &
      penalty_factor,communicator,hpsi,hermiticity_defect,face_balance_defect,ok,message)
    type(s_dg_nodal_common_state), intent(inout) :: state
    integer, intent(in) :: fragment_origins(:,:),fragment_sizes(:,:),global_size(3),communicator
    real(8), intent(in) :: grid_spacing(3),penalty_factor
    complex(8), intent(out) :: hpsi(:,:,:,:,:)
    real(8), intent(out) :: hermiticity_defect,face_balance_defect
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: local_value(:,:,:,:),local_normal(:,:,:,:)
    complex(8), allocatable :: neighbor_value(:,:,:,:),neighbor_normal(:,:,:,:)
    complex(8), allocatable :: face_value(:,:,:,:),face_normal(:,:,:,:)
    complex(8), allocatable :: paired_value(:,:,:,:),paired_normal(:,:,:,:)
    real(8), allocatable :: derivative_weights(:)
    real(8) :: cell_volume,face_weight,local_imag,global_imag,local_balance,global_balance
    complex(8) :: energy
    complex(8), allocatable :: local_projected_h(:,:),global_projected_h(:,:)
    integer :: axis,side,iface,nweight,istate,jstate,allocation_status
    logical :: stage_ok,local_stage_ok
#ifdef USE_MPI
    integer :: requests(24),statuses(MPI_STATUS_SIZE,24),request_count,ierr,opposite
#endif
    ok=state%initialized .and. all(shape(hpsi)==shape(state%psi_core)) .and. &
      all(grid_spacing>0d0) .and. penalty_factor>0d0 .and. &
      all(fragment_sizes(:,state%fragment)==state%core_size)
    call dg_dc_collective_and(ok,communicator,stage_ok)
    if(.not.stage_ok) then
      hpsi=(0d0,0d0); hermiticity_defect=huge(1d0); face_balance_defect=huge(1d0)
      message='DG DC GS adapter: invalid SIPG operator context'
      ok=.false.
      return
    end if
    hpsi=(0d0,0d0)
    cell_volume=product(grid_spacing)
    local_stage_ok=.true.
    do axis=1,3
      nweight=min(5,state%core_size(axis))
      if(nweight<2 .and. global_size(axis)>1) then
        local_stage_ok=.false.
      end if
      if(nweight<2) cycle
      allocate(derivative_weights(nweight))
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        call build_dg_dc_one_sided_derivative_weights(grid_spacing(axis),side,derivative_weights,stage_ok,message)
        local_stage_ok=local_stage_ok .and. stage_ok
        if(.not.stage_ok) cycle
        call pack_face_trace(state,axis,side,derivative_weights,state%faces(iface)%send_value(1,:,:,:,:), &
          state%faces(iface)%send_normal(1,:,:,:,:))
      end do
      deallocate(derivative_weights)
    end do
    call dg_dc_collective_and(local_stage_ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; message='DG DC GS adapter: SIPG trace construction failed collectively'; return
    end if
#ifdef USE_MPI
    request_count=0
    requests=MPI_REQUEST_NULL
    local_stage_ok=.true.
    do axis=1,3
      if(state%core_size(axis)<2) cycle
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        opposite=dg_nodal_face_slot(axis,-side)
        request_count=request_count+1
        call MPI_Irecv(state%faces(iface)%recv_value,size(state%faces(iface)%recv_value),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,100+2*opposite,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
        request_count=request_count+1
        call MPI_Irecv(state%faces(iface)%recv_normal,size(state%faces(iface)%recv_normal),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,101+2*opposite,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
        request_count=request_count+1
        call MPI_Isend(state%faces(iface)%send_value,size(state%faces(iface)%send_value),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,100+2*iface,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
        request_count=request_count+1
        call MPI_Isend(state%faces(iface)%send_normal,size(state%faces(iface)%send_normal),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,101+2*iface,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
      end do
    end do
    call MPI_Waitall(request_count,requests,statuses,ierr)
    local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
    call dg_dc_collective_and(local_stage_ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; message='DG DC GS adapter: SIPG trace exchange failed collectively'; return
    end if
#else
    do iface=1,size(state%faces)
      state%faces(iface)%recv_value=state%faces(iface)%send_value
      state%faces(iface)%recv_normal=state%faces(iface)%send_normal
    end do
#endif
    local_balance=0d0
    local_stage_ok=.true.
    do axis=1,3
      nweight=min(5,state%core_size(axis))
      if(nweight<2) cycle
      allocate(derivative_weights(nweight))
      face_weight=cell_volume/grid_spacing(axis)
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        call build_dg_dc_one_sided_derivative_weights(grid_spacing(axis),side,derivative_weights,stage_ok,message)
        local_stage_ok=local_stage_ok .and. stage_ok
        if(.not.stage_ok) cycle
        allocate(local_value,mold=state%faces(iface)%send_value(1,:,:,:,:))
        allocate(local_normal,mold=local_value)
        allocate(neighbor_value,mold=local_value)
        allocate(neighbor_normal,mold=local_value)
        allocate(face_value,mold=local_value)
        allocate(face_normal,mold=local_value)
        allocate(paired_value,mold=local_value)
        allocate(paired_normal,mold=local_value)
        local_value=state%faces(iface)%send_value(1,:,:,:,:)
        local_normal=state%faces(iface)%send_normal(1,:,:,:,:)
        neighbor_value=state%faces(iface)%recv_value(1,:,:,:,:)
        neighbor_normal=state%faces(iface)%recv_normal(1,:,:,:,:)
        call evaluate_dg_dc_local_sipg_face(local_value,local_normal,neighbor_value,neighbor_normal, &
          grid_spacing(axis),face_weight,penalty_factor,side,face_value,face_normal,stage_ok,message)
        local_stage_ok=local_stage_ok .and. stage_ok
        if(stage_ok) call evaluate_dg_dc_local_sipg_face(neighbor_value,neighbor_normal,local_value,local_normal, &
          grid_spacing(axis),face_weight,penalty_factor,-side,paired_value,paired_normal,stage_ok,message)
        local_stage_ok=local_stage_ok .and. stage_ok
        if(stage_ok) then
          state%faces(iface)%send_value(1,:,:,:,:)=face_value
          state%faces(iface)%send_normal(1,:,:,:,:)=face_normal
        end if
        face_value=face_value/cell_volume
        face_normal=face_normal/cell_volume
        if(stage_ok) call lift_dg_dc_sipg_face(axis,side,derivative_weights,face_value,face_normal,hpsi,stage_ok,message)
        local_stage_ok=local_stage_ok .and. stage_ok
        deallocate(local_value,local_normal,neighbor_value,neighbor_normal,face_value,face_normal, &
          paired_value,paired_normal)
      end do
      deallocate(derivative_weights)
    end do
    call dg_dc_collective_and(local_stage_ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; hpsi=(0d0,0d0); message='DG DC GS adapter: SIPG evaluation/lift failed collectively'; return
    end if
#ifdef USE_MPI
    request_count=0
    requests=MPI_REQUEST_NULL
    local_stage_ok=.true.
    do axis=1,3
      if(state%core_size(axis)<2) cycle
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        opposite=dg_nodal_face_slot(axis,-side)
        request_count=request_count+1
        call MPI_Irecv(state%faces(iface)%recv_value,size(state%faces(iface)%recv_value),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,300+2*opposite,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
        request_count=request_count+1
        call MPI_Irecv(state%faces(iface)%recv_normal,size(state%faces(iface)%recv_normal),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,301+2*opposite,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
        request_count=request_count+1
        call MPI_Isend(state%faces(iface)%send_value,size(state%faces(iface)%send_value),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,300+2*iface,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
        request_count=request_count+1
        call MPI_Isend(state%faces(iface)%send_normal,size(state%faces(iface)%send_normal),MPI_DOUBLE_COMPLEX, &
          state%faces(iface)%neighbor_rank,301+2*iface,communicator,requests(request_count),ierr)
        local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
      end do
    end do
    call MPI_Waitall(request_count,requests,statuses,ierr)
    local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
    call dg_dc_collective_and(local_stage_ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; hpsi=(0d0,0d0); message='DG DC GS adapter: SIPG balance exchange failed collectively'; return
    end if
#else
    do iface=1,size(state%faces)
      state%faces(iface)%recv_value=state%faces(iface)%send_value
      state%faces(iface)%recv_normal=state%faces(iface)%send_normal
    end do
#endif
    local_balance=0d0
    do axis=1,3
      if(state%core_size(axis)<2) cycle
      do side=-1,1,2
        iface=dg_nodal_face_slot(axis,side)
        call measure_dg_dc_face_action_balance(state%faces(iface)%send_value(1,:,:,:,:), &
          state%faces(iface)%send_normal(1,:,:,:,:),state%faces(iface)%recv_value(1,:,:,:,:), &
          state%faces(iface)%recv_normal(1,:,:,:,:), &
          face_weight,stage_ok)
        if(stage_ok) local_balance=max(local_balance,face_weight)
      end do
    end do
    energy=sum(conjg(state%psi_core)*hpsi)
    local_imag=abs(aimag(energy))*cell_volume
    allocate(local_projected_h(state%nstate,state%nstate),global_projected_h(state%nstate,state%nstate), &
      stat=allocation_status)
    local_stage_ok=allocation_status==0
    call dg_dc_collective_and(local_stage_ok,communicator,stage_ok)
    if(stage_ok) then
      do jstate=1,state%nstate
        do istate=1,state%nstate
          local_projected_h(istate,jstate)=sum(conjg(state%psi_core(:,:,:,istate,1))* &
            hpsi(:,:,:,jstate,1))*cell_volume
        end do
      end do
#ifdef USE_MPI
      call MPI_Allreduce(local_projected_h,global_projected_h,size(local_projected_h),MPI_DOUBLE_COMPLEX, &
        MPI_SUM,communicator,ierr)
      local_stage_ok=ierr==MPI_SUCCESS
#else
      global_projected_h=local_projected_h
      local_stage_ok=.true.
#endif
      if(local_stage_ok) then
        face_weight=maxval(abs(global_projected_h-conjg(transpose(global_projected_h))))/ &
          max(1d0,maxval(abs(global_projected_h)))
      else
        face_weight=huge(1d0)
      end if
    else
      face_weight=huge(1d0)
      local_stage_ok=.false.
    end if
    if(allocated(local_projected_h)) deallocate(local_projected_h)
    if(allocated(global_projected_h)) deallocate(global_projected_h)
#ifdef USE_MPI
    call MPI_Allreduce(local_imag,global_imag,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
    call MPI_Allreduce(local_balance,global_balance,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    local_stage_ok=local_stage_ok .and. ierr==MPI_SUCCESS
    call dg_dc_collective_and(local_stage_ok,communicator,ok)
#else
    global_imag=local_imag
    global_balance=local_balance
    ok=.true.
#endif
    ok=ok .and. local_stage_ok
    hermiticity_defect=global_imag
    if(local_stage_ok) hermiticity_defect=max(hermiticity_defect,face_weight)
    face_balance_defect=global_balance
    if(ok) then
      message=''
    else
      message='DG DC GS adapter: SIPG diagnostic reduction failed'
    end if
  end subroutine apply_dg_dc_sipg_operator_mpi

  subroutine pack_face_trace(state,axis,side,derivative_weights,value,normal)
    type(s_dg_nodal_common_state), intent(in) :: state
    integer, intent(in) :: axis,side
    real(8), intent(in) :: derivative_weights(:)
    complex(8), intent(out) :: value(:,:,:,:),normal(:,:,:,:)
    integer :: i,j,k,n,istate,ispin
    value=(0d0,0d0); normal=(0d0,0d0)
    do ispin=1,state%nspin
      do istate=1,state%nstate
        select case(axis)
        case(1)
          do k=1,state%core_size(3); do j=1,state%core_size(2)
            i=merge(1,state%core_size(1),side<0)
            value(j,k,istate,ispin)=state%psi_core(i,j,k,istate,ispin)
            do n=1,size(derivative_weights)
              i=merge(n,state%core_size(1)-n+1,side<0)
              normal(j,k,istate,ispin)=normal(j,k,istate,ispin)+ &
                derivative_weights(n)*state%psi_core(i,j,k,istate,ispin)
            end do
          end do; end do
        case(2)
          do k=1,state%core_size(3); do i=1,state%core_size(1)
            j=merge(1,state%core_size(2),side<0)
            value(i,k,istate,ispin)=state%psi_core(i,j,k,istate,ispin)
            do n=1,size(derivative_weights)
              j=merge(n,state%core_size(2)-n+1,side<0)
              normal(i,k,istate,ispin)=normal(i,k,istate,ispin)+ &
                derivative_weights(n)*state%psi_core(i,j,k,istate,ispin)
            end do
          end do; end do
        case(3)
          do j=1,state%core_size(2); do i=1,state%core_size(1)
            k=merge(1,state%core_size(3),side<0)
            value(i,j,istate,ispin)=state%psi_core(i,j,k,istate,ispin)
            do n=1,size(derivative_weights)
              k=merge(n,state%core_size(3)-n+1,side<0)
              normal(i,j,istate,ispin)=normal(i,j,istate,ispin)+ &
                derivative_weights(n)*state%psi_core(i,j,k,istate,ispin)
            end do
          end do; end do
        end select
      end do
    end do
  end subroutine pack_face_trace

  subroutine build_dg_dc_interior_volume_action(psi,vlocal,coef_lap0,coef_lap,zero_extended_local, &
      interior_volume,ok,message)
    complex(8), intent(in) :: psi(:,:,:,:,:)
    real(8), intent(in) :: vlocal(:,:,:,:),coef_lap0,coef_lap(:,:)
    complex(8), intent(out) :: zero_extended_local(:,:,:,:,:),interior_volume(:,:,:,:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: axis,offset,ix,iy,iz,jx,jy,jz,istate,ispin
    real(8) :: pair_coefficient
    complex(8) :: difference
    ok=size(coef_lap,2)==3 .and. all(shape(vlocal)==[shape(psi(:,:,:,1,1)),size(psi,5)]) .and. &
      all(shape(zero_extended_local)==shape(psi)) .and. all(shape(interior_volume)==shape(psi)) .and. &
      ieee_is_finite(coef_lap0) .and. all(ieee_is_finite(coef_lap)) .and. all(ieee_is_finite(vlocal))
    if(.not.ok) then
      zero_extended_local=(0d0,0d0); interior_volume=(0d0,0d0)
      message='DG DC GS adapter: invalid interior-volume layout'
      return
    end if
    do ispin=1,size(psi,5)
      do istate=1,size(psi,4)
        zero_extended_local(:,:,:,istate,ispin)= &
          (-0.5d0*coef_lap0+vlocal(:,:,:,ispin))*psi(:,:,:,istate,ispin)
        interior_volume(:,:,:,istate,ispin)=vlocal(:,:,:,ispin)*psi(:,:,:,istate,ispin)
      end do
    end do
    do axis=1,3
      do offset=1,size(coef_lap,1)
        pair_coefficient=0.5d0*coef_lap(offset,axis)
        do ispin=1,size(psi,5)
          do istate=1,size(psi,4)
            do iz=1,size(psi,3)
              do iy=1,size(psi,2)
                do ix=1,size(psi,1)
                  jx=ix; jy=iy; jz=iz
                  select case(axis)
                  case(1); jx=ix+offset
                  case(2); jy=iy+offset
                  case(3); jz=iz+offset
                  end select
                  if(jx>size(psi,1) .or. jy>size(psi,2) .or. jz>size(psi,3)) cycle
                  zero_extended_local(ix,iy,iz,istate,ispin)= &
                    zero_extended_local(ix,iy,iz,istate,ispin)- &
                    pair_coefficient*psi(jx,jy,jz,istate,ispin)
                  zero_extended_local(jx,jy,jz,istate,ispin)= &
                    zero_extended_local(jx,jy,jz,istate,ispin)- &
                    pair_coefficient*psi(ix,iy,iz,istate,ispin)
                  difference=psi(ix,iy,iz,istate,ispin)-psi(jx,jy,jz,istate,ispin)
                  interior_volume(ix,iy,iz,istate,ispin)= &
                    interior_volume(ix,iy,iz,istate,ispin)+pair_coefficient*difference
                  interior_volume(jx,jy,jz,istate,ispin)= &
                    interior_volume(jx,jy,jz,istate,ispin)-pair_coefficient*difference
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    ok=all(ieee_is_finite(real(zero_extended_local,8))) .and. &
      all(ieee_is_finite(aimag(zero_extended_local))) .and. &
      all(ieee_is_finite(real(interior_volume,8))) .and. all(ieee_is_finite(aimag(interior_volume)))
    if(ok) then
      message=''
    else
      zero_extended_local=(0d0,0d0); interior_volume=(0d0,0d0)
      message='DG DC GS adapter: nonfinite interior-volume action'
    end if
  end subroutine build_dg_dc_interior_volume_action

  subroutine measure_dg_dc_face_action_balance(local_value,local_normal,peer_value,peer_normal,defect,ok)
    complex(8), intent(in) :: local_value(:,:,:,:),local_normal(:,:,:,:)
    complex(8), intent(in) :: peer_value(:,:,:,:),peer_normal(:,:,:,:)
    real(8), intent(out) :: defect
    logical, intent(out) :: ok
    ok=all(shape(local_value)==shape(local_normal)) .and. all(shape(local_value)==shape(peer_value)) .and. &
      all(shape(local_value)==shape(peer_normal))
    if(ok) then
      defect=max(maxval(abs(local_value+peer_value)),maxval(abs(local_normal+peer_normal)))
      ok=ieee_is_finite(defect)
    else
      defect=huge(1d0)
    end if
  end subroutine measure_dg_dc_face_action_balance

  subroutine measure_dg_dc_cross_hermiticity_mpi(x,hx,y,hy,communicator,defect,ok)
    complex(8), intent(in) :: x(:,:,:),hx(:,:,:),y(:,:,:),hy(:,:,:)
    integer, intent(in) :: communicator
    real(8), intent(out) :: defect
    logical, intent(out) :: ok
    complex(8) :: local_pair(2),global_pair(2)
    ok=all(shape(x)==shape(hx)) .and. all(shape(x)==shape(y)) .and. all(shape(x)==shape(hy))
    call dg_dc_collective_and(ok,communicator,ok)
    if(.not.ok) then
      defect=huge(1d0)
      return
    end if
    local_pair(1)=sum(conjg(x)*hy)
    local_pair(2)=sum(conjg(y)*hx)
#ifdef USE_MPI
    block
      integer :: ierr
      call MPI_Allreduce(local_pair,global_pair,2,MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,ierr)
      ok=ierr==MPI_SUCCESS
    end block
#else
    global_pair=local_pair
#endif
    defect=abs(global_pair(1)-conjg(global_pair(2)))
    ok=ok .and. ieee_is_finite(defect)
  end subroutine measure_dg_dc_cross_hermiticity_mpi

  subroutine validate_dg_dc_candidate_memory(core_size,global_candidate_count,nspin,maximum_bytes, &
      estimated_bytes,ok)
    integer, intent(in) :: core_size(3),global_candidate_count,nspin
    real(8), intent(in) :: maximum_bytes
    real(8), intent(out) :: estimated_bytes
    logical, intent(out) :: ok
    real(8) :: state_elements,matrix_elements
    ok=all(core_size>0) .and. global_candidate_count>0 .and. nspin>0 .and. &
      maximum_bytes>0d0 .and. ieee_is_finite(maximum_bytes)
    if(.not.ok) then
      estimated_bytes=huge(1d0)
      return
    end if
    state_elements=product(real(core_size,8))*real(global_candidate_count,8)*real(nspin,8)
    matrix_elements=real(global_candidate_count,8)**2*real(nspin,8)
    estimated_bytes=16d0*(8d0*state_elements+4d0*matrix_elements)
    ok=ieee_is_finite(estimated_bytes) .and. estimated_bytes<=maximum_bytes
  end subroutine validate_dg_dc_candidate_memory

  subroutine dg_dc_collective_and(local_value,communicator,global_value)
    logical, intent(in) :: local_value
    integer, intent(in) :: communicator
    logical, intent(out) :: global_value
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS) global_value=.false.
#else
    global_value=local_value .and. communicator>=0
#endif
  end subroutine dg_dc_collective_and

end module dg_dc_ground_state_adapter
