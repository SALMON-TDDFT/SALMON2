!
!  Copyright 2019-2026 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!
#include "config.h"
module dg_dc_local_basis_ground_state
#ifdef USE_MPI
  use mpi
#endif
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dg_nodal_sipg, only: s_dg_nodal_sipg_action,evaluate_dg_nodal_sipg_face
  implicit none
  private

  type, public :: s_dg_dc_local_basis_layout
    integer :: fragment_id=0
    integer :: local_basis_count=0
    integer :: global_basis_count=0
    integer :: global_band_count=0
    integer(8) :: geometry_fingerprint=0_8
    integer(8) :: operator_fingerprint=0_8
    integer, allocatable :: basis_offsets(:)
    integer, allocatable :: fragment_ids(:)
  end type s_dg_dc_local_basis_layout
  type, public :: s_dg_dc_local_basis_production_state
    logical :: ready=.false.
    integer :: scf_iterations=0
    integer(8) :: geometry_fingerprint=0_8,operator_fingerprint=0_8
    real(8) :: density_residual=huge(1d0),interface_scale=0d0
    complex(8), allocatable :: coefficient_rows(:,:)
    complex(8), allocatable :: full_fragment_basis(:,:),basis_transform(:,:)
    real(8), allocatable :: eigenvalues(:),occupations(:)
    integer, allocatable :: basis_offsets(:),fragment_ids(:)
    integer :: fragment_id=0,local_basis_count=0,global_basis_count=0,global_band_count=0
    integer :: core_size(3)=0,full_spatial_shape(3)=0
  end type s_dg_dc_local_basis_production_state
  type(s_dg_dc_local_basis_production_state), public, save :: dg_dc_local_basis_state

  public :: initialize_dg_dc_local_basis_layout
  public :: release_dg_dc_local_basis_layout
  public :: assemble_dg_dc_local_basis_sipg_pair
  public :: exchange_dg_dc_local_basis_axis_traces
  public :: pack_dg_dc_local_basis_face_trace
  public :: accumulate_dg_dc_local_basis_interface_rows
  public :: build_dg_dc_six_face_neighbors
  public :: diagnose_dg_dc_six_face_balance
  public :: assemble_dg_dc_local_basis_interface_rows
  public :: compose_dg_dc_local_basis_pair_hamiltonian
  public :: solve_dg_dc_local_basis_bands_reference
  ! Production coefficient CG requires fragment-core-orthonormalized basis rows;
  ! identity_overlap_rows is a distributed identity certificate, not a general S.
  public :: solve_dg_dc_local_basis_bands_cg
  public :: initialize_dg_dc_local_basis_coefficients
  public :: assign_dg_dc_local_basis_occupations
  public :: reconstruct_dg_dc_local_basis_density
  public :: validate_dg_dc_local_basis_density
  public :: diagnose_dg_dc_local_basis_continuation
  public :: orthonormalize_dg_dc_fragment_core_basis
  public :: transform_dg_dc_fragment_buffer_basis
  public :: project_dg_dc_local_basis_volume
  public :: compose_dg_dc_distributed_hamiltonian_rows

contains

  subroutine diagnose_dg_dc_local_basis_continuation(layout,hamiltonian_rows,coefficient_rows, &
      reference_coefficient_rows,occupations,communicator,projector_overlap,orthogonality_defect, &
      hermiticity_defect,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(in) :: hamiltonian_rows(:,:),coefficient_rows(:,:),reference_coefficient_rows(:,:)
    real(8), intent(in) :: occupations(:)
    integer, intent(in) :: communicator
    real(8), intent(out) :: projector_overlap,orthogonality_defect,hermiticity_defect
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: send_buffer(:),receive_buffer(:),local_gram(:,:),global_gram(:,:), &
      local_overlap(:,:),global_overlap(:,:)
    integer, allocatable :: send_counts(:),receive_counts(:),send_displacements(:),receive_displacements(:), &
      occupied_indices(:)
    integer :: nproc,column,noccupied,info,allocation_status,peer,first,last,offset,peer_count
    real(8) :: local_scale,global_scale,local_hermiticity,global_hermiticity

    projector_overlap=0d0;orthogonality_defect=huge(1d0);hermiticity_defect=huge(1d0)
    nproc=size(layout%fragment_ids)
    allocate(send_counts(nproc),receive_counts(nproc),send_displacements(nproc),receive_displacements(nproc), &
      stat=allocation_status)
    ok=allocation_status==0
    if(ok) then
      ok=all(shape(hamiltonian_rows)==[layout%local_basis_count,layout%global_basis_count]) .and. &
        all(shape(coefficient_rows)==[layout%local_basis_count,layout%global_band_count]) .and. &
        all(shape(reference_coefficient_rows)==shape(coefficient_rows)) .and. &
        size(occupations)==layout%global_band_count .and. all(occupations>=0d0)
    end if
    call collective_logical_and(ok,communicator,ok)
    if(.not.ok) then
      message='DG DC local basis: invalid continuation diagnostic layout'; return
    end if
    offset=0
    do peer=1,nproc
      peer_count=layout%basis_offsets(peer)-layout%basis_offsets(peer-1)
      send_counts(peer)=layout%local_basis_count*peer_count
      receive_counts(peer)=send_counts(peer)
      send_displacements(peer)=offset
      receive_displacements(peer)=offset
      offset=offset+send_counts(peer)
    end do
    allocate(send_buffer(offset),receive_buffer(offset), &
      local_gram(layout%global_band_count,layout%global_band_count), &
      global_gram(layout%global_band_count,layout%global_band_count), &
      local_overlap(layout%global_band_count,layout%global_band_count), &
      global_overlap(layout%global_band_count,layout%global_band_count),stat=allocation_status)
    call collective_logical_and(allocation_status==0,communicator,ok)
    if(.not.ok) then
      message='DG DC local basis: continuation diagnostic allocation failed'; return
    end if
    do peer=1,nproc
      first=layout%basis_offsets(peer-1)+1
      last=layout%basis_offsets(peer)
      offset=send_displacements(peer)
      send_buffer(offset+1:offset+send_counts(peer))=reshape(hamiltonian_rows(:,first:last),[send_counts(peer)])
    end do
#ifdef USE_MPI
    call MPI_Alltoallv(send_buffer,send_counts,send_displacements,MPI_DOUBLE_COMPLEX,receive_buffer, &
      receive_counts,receive_displacements,MPI_DOUBLE_COMPLEX,communicator,info)
#else
    receive_buffer=send_buffer
    info=merge(0,1,communicator>=0)
#endif
    call collective_logical_and(info==0,communicator,ok)
    if(.not.ok) then
      message='DG DC local basis: distributed Hermiticity exchange failed'; return
    end if
    local_scale=max(1d0,maxval(abs(hamiltonian_rows)))
    local_hermiticity=0d0
    do peer=1,nproc
      first=layout%basis_offsets(peer-1)+1
      last=layout%basis_offsets(peer)
      peer_count=last-first+1
      offset=receive_displacements(peer)
      local_hermiticity=max(local_hermiticity,maxval(abs(hamiltonian_rows(:,first:last)- &
        conjg(transpose(reshape(receive_buffer(offset+1:offset+receive_counts(peer)), &
        [peer_count,layout%local_basis_count]))))))
    end do
    local_gram=matmul(conjg(transpose(coefficient_rows)),coefficient_rows)
    local_overlap=matmul(conjg(transpose(reference_coefficient_rows)),coefficient_rows)
#ifdef USE_MPI
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,info)
    call collective_logical_and(info==MPI_SUCCESS,communicator,ok)
    if(ok) call MPI_Allreduce(local_hermiticity,global_hermiticity,1,MPI_DOUBLE_PRECISION, &
      MPI_MAX,communicator,info)
    call collective_logical_and(ok .and. info==MPI_SUCCESS,communicator,ok)
    if(ok) call MPI_Allreduce(local_gram,global_gram,size(local_gram),MPI_DOUBLE_COMPLEX, &
      MPI_SUM,communicator,info)
    call collective_logical_and(ok .and. info==MPI_SUCCESS,communicator,ok)
    if(ok) call MPI_Allreduce(local_overlap,global_overlap,size(local_overlap), &
      MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,info)
    call collective_logical_and(ok .and. info==MPI_SUCCESS,communicator,ok)
#else
    global_scale=local_scale
    global_hermiticity=local_hermiticity
    global_gram=local_gram
    global_overlap=local_overlap
#endif
    ok=ok .and. all(ieee_is_finite(real(global_gram,8))) .and. &
      all(ieee_is_finite(aimag(global_gram))) .and. all(ieee_is_finite(real(global_overlap,8))) .and. &
      all(ieee_is_finite(aimag(global_overlap)))
    if(ok) then
      hermiticity_defect=global_hermiticity/global_scale
      do column=1,layout%global_band_count
        global_gram(column,column)=global_gram(column,column)-1d0
      end do
      orthogonality_defect=maxval(abs(global_gram))
      noccupied=count(occupations>0d0)
      ok=noccupied>0
    end if
    if(ok) then
      allocate(occupied_indices(noccupied),stat=allocation_status)
      ok=allocation_status==0
    end if
    if(ok) then
      occupied_indices=pack([(column,column=1,layout%global_band_count)],occupations>0d0)
      projector_overlap=sum(abs(global_overlap(occupied_indices,occupied_indices))**2)/real(noccupied,8)
      ok=all(ieee_is_finite([projector_overlap,orthogonality_defect,hermiticity_defect])) .and. &
        projector_overlap>=0d0 .and. projector_overlap<=1d0+1d-10
    end if
    call collective_logical_and(ok,communicator,ok)
    if(ok) then
      message=''
    else
      message='DG DC local basis: invalid continuation diagnostics'
    end if
  end subroutine diagnose_dg_dc_local_basis_continuation

  subroutine transform_dg_dc_fragment_buffer_basis(raw_basis,basis_transform,effective_basis_count, &
      transformed_basis,ok,message)
    complex(8), intent(in) :: raw_basis(:,:),basis_transform(:,:)
    integer, intent(in) :: effective_basis_count
    complex(8), intent(out) :: transformed_basis(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message

    transformed_basis=(0d0,0d0)
    ok=size(raw_basis,1)>0 .and. size(raw_basis,2)>0 .and. &
      all(shape(basis_transform)==[size(raw_basis,2),size(raw_basis,2)]) .and. &
      effective_basis_count>0 .and. effective_basis_count<=size(raw_basis,2) .and. &
      all(shape(transformed_basis)==[size(raw_basis,1),effective_basis_count]) .and. &
      all(ieee_is_finite(real(raw_basis,8))) .and. all(ieee_is_finite(aimag(raw_basis))) .and. &
      all(ieee_is_finite(real(basis_transform,8))) .and. all(ieee_is_finite(aimag(basis_transform)))
    if(ok) then
      transformed_basis=matmul(raw_basis,basis_transform(:,1:effective_basis_count))
      ok=all(ieee_is_finite(real(transformed_basis,8))) .and. &
        all(ieee_is_finite(aimag(transformed_basis)))
    end if
    if(ok) then
      message=''
    else
      transformed_basis=(0d0,0d0)
      message='DG DC local basis: invalid full fragment basis transform'
    end if
  end subroutine transform_dg_dc_fragment_buffer_basis

  subroutine project_dg_dc_local_basis_volume(basis,hbasis,quadrature_weight,volume_block,ok,message)
    complex(8), intent(in) :: basis(:,:),hbasis(:,:)
    real(8), intent(in) :: quadrature_weight
    complex(8), intent(out) :: volume_block(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message

    volume_block=(0d0,0d0)
    ok=size(basis,1)>0 .and. size(basis,2)>0 .and. all(shape(hbasis)==shape(basis)) .and. &
      all(shape(volume_block)==[size(basis,2),size(basis,2)]) .and. &
      quadrature_weight>0d0 .and. ieee_is_finite(quadrature_weight) .and. &
      all(ieee_is_finite(real(basis,8))) .and. all(ieee_is_finite(aimag(basis))) .and. &
      all(ieee_is_finite(real(hbasis,8))) .and. all(ieee_is_finite(aimag(hbasis)))
    if(.not.ok) then
      message='DG DC local basis: invalid DC volume projection'
      return
    end if
    volume_block=quadrature_weight*matmul(conjg(transpose(basis)),hbasis)
    ok=all(ieee_is_finite(real(volume_block,8))) .and. all(ieee_is_finite(aimag(volume_block)))
    if(ok) then
      message=''
    else
      volume_block=(0d0,0d0)
      message='DG DC local basis: nonfinite DC volume block'
    end if
  end subroutine project_dg_dc_local_basis_volume

  subroutine initialize_dg_dc_local_basis_coefficients(layout,coefficient_rows,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(out) :: coefficient_rows(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: rank_index,local_row,global_row

    coefficient_rows=(0d0,0d0)
    ok=allocated(layout%basis_offsets) .and. allocated(layout%fragment_ids)
    if(ok) ok=size(layout%basis_offsets)==size(layout%fragment_ids)+1 .and. &
      all(shape(coefficient_rows)==[layout%local_basis_count,layout%global_band_count]) .and. &
      layout%global_basis_count>=layout%global_band_count
    rank_index=-1
    if(ok) then
      do local_row=1,size(layout%fragment_ids)
        if(layout%fragment_ids(local_row)==layout%fragment_id) rank_index=local_row-1
      end do
      ok=rank_index>=0
    end if
    if(ok) ok=layout%basis_offsets(rank_index+1)-layout%basis_offsets(rank_index)== &
      layout%local_basis_count
    if(.not.ok) then
      message='DG DC local basis: invalid coefficient initialization layout'
      return
    end if
    do local_row=1,layout%local_basis_count
      global_row=layout%basis_offsets(rank_index)+local_row
      if(global_row<=layout%global_band_count) coefficient_rows(local_row,global_row)=1d0
    end do
    message=''
  end subroutine initialize_dg_dc_local_basis_coefficients

  subroutine assign_dg_dc_local_basis_occupations(electron_count,occupations,ok,message)
    real(8), intent(in) :: electron_count
    real(8), intent(out) :: occupations(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: fully_occupied
    real(8) :: remainder

    occupations=0d0
    ok=size(occupations)>0 .and. electron_count>=0d0 .and. ieee_is_finite(electron_count) .and. &
      electron_count<=2d0*real(size(occupations),8)+10d0*epsilon(1d0)
    if(.not.ok) then
      message='DG DC local basis: insufficient global bands for electron count'
      return
    end if
    fully_occupied=min(size(occupations),int(electron_count/2d0))
    if(fully_occupied>0) occupations(1:fully_occupied)=2d0
    remainder=electron_count-2d0*real(fully_occupied,8)
    if(remainder>10d0*epsilon(1d0) .and. fully_occupied<size(occupations)) &
      occupations(fully_occupied+1)=remainder
    ok=all(occupations>=0d0) .and. all(occupations<=2d0) .and. &
      abs(sum(occupations)-electron_count)<=10d0*epsilon(1d0)*max(1d0,electron_count)
    if(ok) then
      message=''
    else
      occupations=0d0
      message='DG DC local basis: invalid global occupations'
    end if
  end subroutine assign_dg_dc_local_basis_occupations

  subroutine exchange_dg_dc_local_basis_axis_traces(minus_value,minus_normal,plus_value,plus_normal, &
      minus_rank,plus_rank,tag,communicator,neighbor_minus_value,neighbor_minus_normal, &
      neighbor_plus_value,neighbor_plus_normal,ok,message)
    complex(8), intent(in) :: minus_value(:,:),minus_normal(:,:),plus_value(:,:),plus_normal(:,:)
    integer, intent(in) :: minus_rank,plus_rank,tag,communicator
    complex(8), intent(out) :: neighbor_minus_value(:,:),neighbor_minus_normal(:,:)
    complex(8), intent(out) :: neighbor_plus_value(:,:),neighbor_plus_normal(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: local_shape(2),remote_minus_shape(2),remote_plus_shape(2)
    logical :: local_ok,stage_ok
#ifdef USE_MPI
    integer :: ierr,nproc
#endif

    neighbor_minus_value=(0d0,0d0); neighbor_minus_normal=(0d0,0d0)
    neighbor_plus_value=(0d0,0d0); neighbor_plus_normal=(0d0,0d0)
    local_shape=shape(minus_value)
    local_ok=all(shape(minus_normal)==local_shape) .and. all(shape(plus_value)==local_shape) .and. &
      all(shape(plus_normal)==local_shape) .and. all(local_shape>0) .and. tag>=0 .and. &
      all(ieee_is_finite(real(minus_value,8))) .and. all(ieee_is_finite(aimag(minus_value))) .and. &
      all(ieee_is_finite(real(minus_normal,8))) .and. all(ieee_is_finite(aimag(minus_normal))) .and. &
      all(ieee_is_finite(real(plus_value,8))) .and. all(ieee_is_finite(aimag(plus_value))) .and. &
      all(ieee_is_finite(real(plus_normal,8))) .and. all(ieee_is_finite(aimag(plus_normal)))
#ifdef USE_MPI
    call MPI_Comm_size(communicator,nproc,ierr)
    local_ok=local_ok .and. ierr==MPI_SUCCESS .and. minus_rank>=0 .and. minus_rank<nproc .and. &
      plus_rank>=0 .and. plus_rank<nproc
    call collective_logical_and(local_ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; message='DG DC local basis: invalid axis trace exchange metadata'; return
    end if
    call MPI_Sendrecv(local_shape,2,MPI_INTEGER,plus_rank,tag,remote_minus_shape,2,MPI_INTEGER, &
      minus_rank,tag,communicator,MPI_STATUS_IGNORE,ierr)
    local_ok=ierr==MPI_SUCCESS
    call MPI_Sendrecv(local_shape,2,MPI_INTEGER,minus_rank,tag+1,remote_plus_shape,2,MPI_INTEGER, &
      plus_rank,tag+1,communicator,MPI_STATUS_IGNORE,ierr)
    local_ok=local_ok .and. ierr==MPI_SUCCESS .and. &
      all(shape(neighbor_minus_value)==remote_minus_shape) .and. &
      all(shape(neighbor_minus_normal)==remote_minus_shape) .and. &
      all(shape(neighbor_plus_value)==remote_plus_shape) .and. &
      all(shape(neighbor_plus_normal)==remote_plus_shape) .and. &
      remote_minus_shape(1)==local_shape(1) .and. remote_plus_shape(1)==local_shape(1)
    call collective_logical_and(local_ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; message='DG DC local basis: inconsistent neighbor axis trace shape'; return
    end if
    call MPI_Sendrecv(plus_value,size(plus_value),MPI_DOUBLE_COMPLEX,plus_rank,tag+2, &
      neighbor_minus_value,size(neighbor_minus_value),MPI_DOUBLE_COMPLEX,minus_rank,tag+2, &
      communicator,MPI_STATUS_IGNORE,ierr)
    local_ok=ierr==MPI_SUCCESS
    call MPI_Sendrecv(minus_value,size(minus_value),MPI_DOUBLE_COMPLEX,minus_rank,tag+3, &
      neighbor_plus_value,size(neighbor_plus_value),MPI_DOUBLE_COMPLEX,plus_rank,tag+3, &
      communicator,MPI_STATUS_IGNORE,ierr)
    local_ok=local_ok .and. ierr==MPI_SUCCESS
    call MPI_Sendrecv(plus_normal,size(plus_normal),MPI_DOUBLE_COMPLEX,plus_rank,tag+4, &
      neighbor_minus_normal,size(neighbor_minus_normal),MPI_DOUBLE_COMPLEX,minus_rank,tag+4, &
      communicator,MPI_STATUS_IGNORE,ierr)
    local_ok=local_ok .and. ierr==MPI_SUCCESS
    call MPI_Sendrecv(minus_normal,size(minus_normal),MPI_DOUBLE_COMPLEX,minus_rank,tag+5, &
      neighbor_plus_normal,size(neighbor_plus_normal),MPI_DOUBLE_COMPLEX,plus_rank,tag+5, &
      communicator,MPI_STATUS_IGNORE,ierr)
    local_ok=local_ok .and. ierr==MPI_SUCCESS
    call collective_logical_and(local_ok,communicator,ok)
#else
    ok=communicator>=0 .and. minus_rank==0 .and. plus_rank==0 .and. local_ok .and. &
      all(shape(neighbor_minus_value)==local_shape) .and. all(shape(neighbor_minus_normal)==local_shape) .and. &
      all(shape(neighbor_plus_value)==local_shape) .and. all(shape(neighbor_plus_normal)==local_shape)
    if(ok) then
      neighbor_minus_value=plus_value; neighbor_minus_normal=plus_normal
      neighbor_plus_value=minus_value; neighbor_plus_normal=minus_normal
    end if
#endif
    if(ok) then
      ok=all(ieee_is_finite(real(neighbor_minus_value,8))) .and. &
        all(ieee_is_finite(aimag(neighbor_minus_value))) .and. &
        all(ieee_is_finite(real(neighbor_minus_normal,8))) .and. &
        all(ieee_is_finite(aimag(neighbor_minus_normal))) .and. &
        all(ieee_is_finite(real(neighbor_plus_value,8))) .and. &
        all(ieee_is_finite(aimag(neighbor_plus_value))) .and. &
        all(ieee_is_finite(real(neighbor_plus_normal,8))) .and. &
        all(ieee_is_finite(aimag(neighbor_plus_normal)))
    end if
    if(ok) then
      message=''
    else
      neighbor_minus_value=(0d0,0d0); neighbor_minus_normal=(0d0,0d0)
      neighbor_plus_value=(0d0,0d0); neighbor_plus_normal=(0d0,0d0)
      message='DG DC local basis: axis trace exchange failed collectively'
    end if
  end subroutine exchange_dg_dc_local_basis_axis_traces

  subroutine pack_dg_dc_local_basis_face_trace(basis,axis,side,derivative_weights,value,normal,ok,message)
    complex(8), intent(in) :: basis(:,:,:,:)
    integer, intent(in) :: axis,side
    real(8), intent(in) :: derivative_weights(:)
    complex(8), intent(out) :: value(:,:),normal(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: ngrid(3),nface,nbasis,point,ibasis,n,ix,iy,iz

    value=(0d0,0d0); normal=(0d0,0d0)
    ngrid=[size(basis,1),size(basis,2),size(basis,3)]
    nbasis=size(basis,4)
    nface=product(ngrid)/ngrid(max(1,min(3,axis)))
    ok=axis>=1 .and. axis<=3 .and. abs(side)==1 .and. nbasis>0 .and. &
      size(derivative_weights)>=2 .and. size(derivative_weights)<=ngrid(max(1,min(3,axis))) .and. &
      all(shape(value)==[nface,nbasis]) .and. all(shape(normal)==shape(value)) .and. &
      all(ieee_is_finite(real(basis,8))) .and. all(ieee_is_finite(aimag(basis))) .and. &
      all(ieee_is_finite(derivative_weights))
    if(.not.ok) then
      message='DG DC local basis: invalid face trace packing layout'
      return
    end if
    do ibasis=1,nbasis
      point=0
      select case(axis)
      case(1)
        do iz=1,ngrid(3); do iy=1,ngrid(2)
          point=point+1
          ix=merge(1,ngrid(1),side<0)
          value(point,ibasis)=basis(ix,iy,iz,ibasis)
          do n=1,size(derivative_weights)
            ix=merge(n,ngrid(1)-n+1,side<0)
            normal(point,ibasis)=normal(point,ibasis)+derivative_weights(n)*basis(ix,iy,iz,ibasis)
          end do
        end do; end do
      case(2)
        do iz=1,ngrid(3); do ix=1,ngrid(1)
          point=point+1
          iy=merge(1,ngrid(2),side<0)
          value(point,ibasis)=basis(ix,iy,iz,ibasis)
          do n=1,size(derivative_weights)
            iy=merge(n,ngrid(2)-n+1,side<0)
            normal(point,ibasis)=normal(point,ibasis)+derivative_weights(n)*basis(ix,iy,iz,ibasis)
          end do
        end do; end do
      case(3)
        do iy=1,ngrid(2); do ix=1,ngrid(1)
          point=point+1
          iz=merge(1,ngrid(3),side<0)
          value(point,ibasis)=basis(ix,iy,iz,ibasis)
          do n=1,size(derivative_weights)
            iz=merge(n,ngrid(3)-n+1,side<0)
            normal(point,ibasis)=normal(point,ibasis)+derivative_weights(n)*basis(ix,iy,iz,ibasis)
          end do
        end do; end do
      end select
    end do
    ok=all(ieee_is_finite(real(value,8))) .and. all(ieee_is_finite(aimag(value))) .and. &
      all(ieee_is_finite(real(normal,8))) .and. all(ieee_is_finite(aimag(normal)))
    if(ok) then
      message=''
    else
      value=(0d0,0d0); normal=(0d0,0d0)
      message='DG DC local basis: nonfinite packed face trace'
    end if
  end subroutine pack_dg_dc_local_basis_face_trace

  subroutine accumulate_dg_dc_local_basis_interface_rows(layout,neighbor_fragment,local_block, &
      neighbor_block,interface_rows,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    integer, intent(in) :: neighbor_fragment
    complex(8), intent(in) :: local_block(:,:),neighbor_block(:,:)
    complex(8), intent(inout) :: interface_rows(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: local_rank,neighbor_rank,index,local_first,local_last,neighbor_first,neighbor_last

    local_rank=-1; neighbor_rank=-1
    ok=allocated(layout%basis_offsets) .and. allocated(layout%fragment_ids)
    if(ok) then
      do index=1,size(layout%fragment_ids)
        if(layout%fragment_ids(index)==layout%fragment_id) local_rank=index-1
        if(layout%fragment_ids(index)==neighbor_fragment) neighbor_rank=index-1
      end do
      ok=local_rank>=0 .and. neighbor_rank>=0 .and. &
        size(layout%basis_offsets)==size(layout%fragment_ids)+1
    end if
    if(ok) then
      local_first=layout%basis_offsets(local_rank)+1
      local_last=layout%basis_offsets(local_rank+1)
      neighbor_first=layout%basis_offsets(neighbor_rank)+1
      neighbor_last=layout%basis_offsets(neighbor_rank+1)
      ok=local_last-local_first+1==layout%local_basis_count .and. &
        all(shape(local_block)==[layout%local_basis_count,layout%local_basis_count]) .and. &
        all(shape(neighbor_block)==[layout%local_basis_count,neighbor_last-neighbor_first+1]) .and. &
        all(shape(interface_rows)==[layout%local_basis_count,layout%global_basis_count])
    end if
    if(ok) ok=all(ieee_is_finite(real(local_block,8))) .and. all(ieee_is_finite(aimag(local_block))) .and. &
      all(ieee_is_finite(real(neighbor_block,8))) .and. all(ieee_is_finite(aimag(neighbor_block))) .and. &
      all(ieee_is_finite(real(interface_rows,8))) .and. all(ieee_is_finite(aimag(interface_rows)))
    if(.not.ok) then
      message='DG DC local basis: invalid interface row block'
      return
    end if
    interface_rows(:,local_first:local_last)=interface_rows(:,local_first:local_last)+local_block
    interface_rows(:,neighbor_first:neighbor_last)= &
      interface_rows(:,neighbor_first:neighbor_last)+neighbor_block
    ok=all(ieee_is_finite(real(interface_rows,8))) .and. all(ieee_is_finite(aimag(interface_rows)))
    if(ok) then
      message=''
    else
      message='DG DC local basis: nonfinite accumulated interface rows'
    end if
  end subroutine accumulate_dg_dc_local_basis_interface_rows

  subroutine build_dg_dc_six_face_neighbors(fragment,fragment_origins,fragment_sizes,global_size, &
      neighbors,periodic_shifts,ok,message)
    integer, intent(in) :: fragment,fragment_origins(:,:),fragment_sizes(:,:),global_size(3)
    integer, intent(out) :: neighbors(3,2),periodic_shifts(3,2)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: axis,side_index,side,candidate,tangent,target,match_count
    logical :: matches

    neighbors=0; periodic_shifts=0
    ok=size(fragment_origins,1)==3 .and. all(shape(fragment_sizes)==shape(fragment_origins)) .and. &
      fragment>=1 .and. fragment<=size(fragment_origins,2) .and. all(global_size>0) .and. &
      all(fragment_sizes>0) .and. all(fragment_origins>=0)
    if(ok) ok=all(fragment_origins+fragment_sizes<=spread(global_size,2,size(fragment_origins,2)))
    if(.not.ok) then
      message='DG DC local basis: invalid six-face topology'
      return
    end if
    do axis=1,3
      do side_index=1,2
        side=merge(-1,1,side_index==1)
        if(side<0) then
          target=fragment_origins(axis,fragment)
          if(target==0) periodic_shifts(axis,side_index)=-1
        else
          target=modulo(fragment_origins(axis,fragment)+fragment_sizes(axis,fragment),global_size(axis))
          if(fragment_origins(axis,fragment)+fragment_sizes(axis,fragment)==global_size(axis)) &
            periodic_shifts(axis,side_index)=1
        end if
        match_count=0
        do candidate=1,size(fragment_origins,2)
          matches=.true.
          do tangent=1,3
            if(tangent==axis) cycle
            matches=matches .and. fragment_origins(tangent,candidate)==fragment_origins(tangent,fragment) .and. &
              fragment_sizes(tangent,candidate)==fragment_sizes(tangent,fragment)
          end do
          if(side<0) then
            matches=matches .and. modulo(fragment_origins(axis,candidate)+fragment_sizes(axis,candidate), &
              global_size(axis))==target
          else
            matches=matches .and. fragment_origins(axis,candidate)==target
          end if
          if(matches) then
            match_count=match_count+1
            neighbors(axis,side_index)=candidate
          end if
        end do
        if(match_count/=1) ok=.false.
      end do
    end do
    if(ok) then
      message=''
    else
      neighbors=0; periodic_shifts=0
      message='DG DC local basis: missing or ambiguous six-face neighbor'
    end if
  end subroutine build_dg_dc_six_face_neighbors

  subroutine diagnose_dg_dc_six_face_balance(fragment_origins,fragment_sizes,global_size,communicator, &
      defect,ok,message)
    integer, intent(in) :: fragment_origins(:,:),fragment_sizes(:,:),global_size(3),communicator
    real(8), intent(out) :: defect
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: fragment,axis,side_index,opposite,neighbor,canonical_count,expected_count
    integer :: neighbors(3,2),shifts(3,2),peer_neighbors(3,2),peer_shifts(3,2)
    logical :: local_ok,stage_ok
    real(8) :: local_defect

    local_ok=size(fragment_origins,1)==3 .and. all(shape(fragment_sizes)==shape(fragment_origins)) .and. &
      size(fragment_origins,2)>0
    local_defect=0d0
    canonical_count=0
    if(local_ok) then
      do fragment=1,size(fragment_origins,2)
        call build_dg_dc_six_face_neighbors(fragment,fragment_origins,fragment_sizes,global_size, &
          neighbors,shifts,stage_ok,message)
        if(.not.stage_ok) then; local_ok=.false.; exit; end if
        do axis=1,3
          do side_index=1,2
            neighbor=neighbors(axis,side_index)
            opposite=3-side_index
            call build_dg_dc_six_face_neighbors(neighbor,fragment_origins,fragment_sizes,global_size, &
              peer_neighbors,peer_shifts,stage_ok,message)
            if(.not.stage_ok) then; local_ok=.false.; exit; end if
            if(peer_neighbors(axis,opposite)/=fragment .or. &
              peer_shifts(axis,opposite)/=-shifts(axis,side_index)) local_defect=1d0
            if(fragment<neighbor .or. (fragment==neighbor .and. side_index==2)) &
              canonical_count=canonical_count+1
          end do
          if(.not.local_ok) exit
        end do
        if(.not.local_ok) exit
      end do
    end if
    expected_count=3*size(fragment_origins,2)
    if(canonical_count/=expected_count) local_defect=1d0
    local_ok=local_ok .and. local_defect==0d0
    call collective_logical_and(local_ok,communicator,ok)
#ifdef USE_MPI
    block
      integer :: ierr
      call MPI_Allreduce(local_defect,defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
      ok=ok .and. ierr==MPI_SUCCESS
    end block
#else
    defect=local_defect
#endif
    if(ok) then
      message=''
    else
      defect=huge(1d0)
      message='DG DC local basis: unbalanced canonical six-face schedule'
    end if
  end subroutine diagnose_dg_dc_six_face_balance

  subroutine assemble_dg_dc_local_basis_interface_rows(layout,basis,fragment_origins,fragment_sizes, &
      global_size,grid_spacing,penalty_factor,communicator,interface_rows,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(in) :: basis(:,:,:,:)
    integer, intent(in) :: fragment_origins(:,:),fragment_sizes(:,:),global_size(3),communicator
    real(8), intent(in) :: grid_spacing(3),penalty_factor
    complex(8), intent(out) :: interface_rows(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: minus_value(:,:),minus_normal(:,:),plus_value(:,:),plus_normal(:,:)
    complex(8), allocatable :: neighbor_minus_value(:,:),neighbor_minus_normal(:,:)
    complex(8), allocatable :: neighbor_plus_value(:,:),neighbor_plus_normal(:,:)
    complex(8), allocatable :: hll(:,:),hln(:,:),hnl(:,:),hnn(:,:)
    real(8), allocatable :: derivative_weights(:)
    integer :: neighbors(3,2),periodic_shifts(3,2),shift_vector(3)
    integer :: axis,side_index,nweight,nface,neighbor,neighbor_rank_index,neighbor_count,index
    real(8) :: face_weight
    logical :: stage_ok,all_stage_ok

    interface_rows=(0d0,0d0)
    stage_ok=layout%fragment_id>=1 .and. layout%fragment_id<=size(fragment_origins,2) .and. &
      layout%local_basis_count==size(basis,4) .and. &
      all(shape(interface_rows)==[layout%local_basis_count,layout%global_basis_count]) .and. &
      all(shape(fragment_origins)==shape(fragment_sizes)) .and. size(fragment_origins,1)==3 .and. &
      all(grid_spacing>0d0) .and. penalty_factor>0d0 .and. all(ieee_is_finite(grid_spacing)) .and. &
      ieee_is_finite(penalty_factor) .and. all(ieee_is_finite(real(basis,8))) .and. &
      all(ieee_is_finite(aimag(basis)))
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(.not.all_stage_ok) then
      ok=.false.; message='DG DC local basis: invalid interface assembly context'; return
    end if
    call build_dg_dc_six_face_neighbors(layout%fragment_id,fragment_origins,fragment_sizes,global_size, &
      neighbors,periodic_shifts,stage_ok,message)
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(.not.all_stage_ok) then
      ok=.false.; message='DG DC local basis: six-face topology failed collectively'; return
    end if
    do axis=1,3
      if(all(neighbors(axis,:)==layout%fragment_id)) cycle
      nweight=min(5,size(basis,axis))
      stage_ok=nweight>=2
      call collective_logical_and(stage_ok,communicator,all_stage_ok)
      if(.not.all_stage_ok) then
        ok=.false.; message='DG DC local basis: insufficient one-sided derivative points'; return
      end if
      nface=product([size(basis,1),size(basis,2),size(basis,3)])/size(basis,axis)
      allocate(derivative_weights(nweight))
      call build_dg_dc_local_basis_outward_weights(grid_spacing(axis),derivative_weights,stage_ok)
      allocate(minus_value(nface,layout%local_basis_count),minus_normal(nface,layout%local_basis_count), &
        plus_value(nface,layout%local_basis_count),plus_normal(nface,layout%local_basis_count))
      call pack_dg_dc_local_basis_face_trace(basis,axis,-1,derivative_weights,minus_value,minus_normal, &
        stage_ok,message)
      if(stage_ok) call pack_dg_dc_local_basis_face_trace(basis,axis,1,derivative_weights,plus_value,plus_normal, &
        stage_ok,message)
      call collective_logical_and(stage_ok,communicator,all_stage_ok)
      if(.not.all_stage_ok) then
        ok=.false.; message='DG DC local basis: face packing failed collectively'; return
      end if
      neighbor_rank_index=-1
      do index=1,size(layout%fragment_ids)
        if(layout%fragment_ids(index)==neighbors(axis,1)) neighbor_rank_index=index-1
      end do
      stage_ok=neighbor_rank_index>=0
      if(stage_ok) neighbor_count=layout%basis_offsets(neighbor_rank_index+1)- &
        layout%basis_offsets(neighbor_rank_index)
      call collective_logical_and(stage_ok,communicator,all_stage_ok)
      if(.not.all_stage_ok) then
        ok=.false.; message='DG DC local basis: minus-face neighbor layout failed'; return
      end if
      allocate(neighbor_minus_value(nface,neighbor_count),neighbor_minus_normal(nface,neighbor_count))
      neighbor_rank_index=-1
      do index=1,size(layout%fragment_ids)
        if(layout%fragment_ids(index)==neighbors(axis,2)) neighbor_rank_index=index-1
      end do
      stage_ok=neighbor_rank_index>=0
      if(stage_ok) neighbor_count=layout%basis_offsets(neighbor_rank_index+1)- &
        layout%basis_offsets(neighbor_rank_index)
      call collective_logical_and(stage_ok,communicator,all_stage_ok)
      if(.not.all_stage_ok) then
        ok=.false.; message='DG DC local basis: plus-face neighbor layout failed'; return
      end if
      allocate(neighbor_plus_value(nface,neighbor_count),neighbor_plus_normal(nface,neighbor_count))
      call exchange_dg_dc_local_basis_axis_traces(minus_value,minus_normal,plus_value,plus_normal, &
        neighbors(axis,1)-1,neighbors(axis,2)-1,1000+10*axis,communicator, &
        neighbor_minus_value,neighbor_minus_normal,neighbor_plus_value,neighbor_plus_normal,stage_ok,message)
      if(.not.stage_ok) then
        ok=.false.; return
      end if
      face_weight=product(grid_spacing)/grid_spacing(axis)
      do side_index=1,2
        neighbor=neighbors(axis,side_index)
        if(side_index==1) then
          neighbor_count=size(neighbor_minus_value,2)
        else
          neighbor_count=size(neighbor_plus_value,2)
        end if
        allocate(hll(layout%local_basis_count,layout%local_basis_count), &
          hln(layout%local_basis_count,neighbor_count),hnl(neighbor_count,layout%local_basis_count), &
          hnn(neighbor_count,neighbor_count))
        shift_vector=0
        shift_vector(axis)=periodic_shifts(axis,side_index)
        if(side_index==1) then
          call assemble_dg_dc_local_basis_sipg_pair(minus_value,minus_normal,neighbor_minus_value, &
            neighbor_minus_normal,grid_spacing(axis),face_weight,penalty_factor,.true.,.true.,.false., &
            shift_vector,shift_vector,neighbor,neighbor,hll,hln,hnl,hnn,stage_ok,message)
        else
          call assemble_dg_dc_local_basis_sipg_pair(plus_value,plus_normal,neighbor_plus_value, &
            neighbor_plus_normal,grid_spacing(axis),face_weight,penalty_factor,.true.,.true.,.false., &
            shift_vector,shift_vector,neighbor,neighbor,hll,hln,hnl,hnn,stage_ok,message)
        end if
        if(stage_ok) call accumulate_dg_dc_local_basis_interface_rows(layout,neighbor,hll,hln, &
          interface_rows,stage_ok,message)
        call collective_logical_and(stage_ok,communicator,all_stage_ok)
        deallocate(hll,hln,hnl,hnn)
        if(.not.all_stage_ok) then
          ok=.false.; message='DG DC local basis: SIPG row assembly failed collectively'; return
        end if
      end do
      deallocate(derivative_weights,minus_value,minus_normal,plus_value,plus_normal, &
        neighbor_minus_value,neighbor_minus_normal,neighbor_plus_value,neighbor_plus_normal)
    end do
    call validate_dg_dc_local_basis_hermitian_rows(layout,interface_rows,communicator,ok)
    if(ok) then
      message=''
    else
      interface_rows=(0d0,0d0)
      message='DG DC local basis: interface rows are not Hermitian'
    end if
  end subroutine assemble_dg_dc_local_basis_interface_rows

  subroutine build_dg_dc_local_basis_outward_weights(grid_spacing,weights,ok)
    real(8), intent(in) :: grid_spacing
    real(8), intent(out) :: weights(:)
    logical, intent(out) :: ok
    weights=0d0
    ok=grid_spacing>0d0 .and. ieee_is_finite(grid_spacing)
    if(.not.ok) return
    select case(size(weights))
    case(2)
      weights=[1d0,-1d0]/grid_spacing
    case(3)
      weights=[3d0,-4d0,1d0]/(2d0*grid_spacing)
    case(4)
      weights=[11d0,-18d0,9d0,-2d0]/(6d0*grid_spacing)
    case(5)
      weights=[25d0,-48d0,36d0,-16d0,3d0]/(12d0*grid_spacing)
    case default
      ok=.false.
    end select
    ok=ok .and. all(ieee_is_finite(weights))
  end subroutine build_dg_dc_local_basis_outward_weights

  subroutine validate_dg_dc_local_basis_hermitian_rows(layout,rows,communicator,ok)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(in) :: rows(:,:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    real(8) :: local_defect,global_defect,local_scale,global_scale
    complex(8), allocatable :: reference_row(:)
    integer :: rank,nproc,global_row,local_row,owner,index,ierr
    logical :: stage_ok

    rank=0; nproc=1
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,rank,ierr)
    ok=ierr==MPI_SUCCESS
    call MPI_Comm_size(communicator,nproc,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
#else
    ok=communicator>=0
#endif
    ok=ok .and. allocated(layout%basis_offsets) .and. &
      all(shape(rows)==[layout%local_basis_count,layout%global_basis_count])
    call collective_logical_and(ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; return
    end if
    local_defect=0d0
    local_scale=max(1d0,maxval(abs(rows)))
    allocate(reference_row(layout%global_basis_count))
    do global_row=1,layout%global_basis_count
      owner=0
      do index=0,nproc-1
        if(global_row>layout%basis_offsets(index) .and. &
          global_row<=layout%basis_offsets(index+1)) owner=index
      end do
      reference_row=(0d0,0d0)
      if(rank==owner) then
        local_row=global_row-layout%basis_offsets(rank)
        reference_row=rows(local_row,:)
      end if
#ifdef USE_MPI
      call MPI_Bcast(reference_row,layout%global_basis_count,MPI_DOUBLE_COMPLEX,owner,communicator,ierr)
      ok=ierr==MPI_SUCCESS
#endif
      if(ok) then
        do local_row=1,layout%local_basis_count
          local_defect=max(local_defect,abs(rows(local_row,global_row)- &
            conjg(reference_row(layout%basis_offsets(rank)+local_row))))
        end do
      end if
    end do
    deallocate(reference_row)
#ifdef USE_MPI
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
#else
    global_defect=local_defect
    global_scale=local_scale
#endif
    ok=ok .and. global_defect<=1d-11*global_scale
    call collective_logical_and(ok,communicator,stage_ok)
    ok=stage_ok
  end subroutine validate_dg_dc_local_basis_hermitian_rows

  subroutine orthonormalize_dg_dc_fragment_core_basis(raw_basis,quadrature_weight,relative_cutoff, &
      orthonormal_basis,basis_transform,effective_basis_count,ok,message)
    complex(8), intent(in) :: raw_basis(:,:)
    real(8), intent(in) :: quadrature_weight,relative_cutoff
    complex(8), intent(out) :: orthonormal_basis(:,:),basis_transform(:,:)
    integer, intent(out) :: effective_basis_count
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: gram(:,:),work(:)
    real(8), allocatable :: eigenvalues(:),rwork(:)
    integer :: nbasis,lwork,info,allocation_status,source_column,output_column
    real(8) :: maximum_eigenvalue,maximum_basis_magnitude,safe_basis_magnitude

    orthonormal_basis=(0d0,0d0)
    basis_transform=(0d0,0d0)
    effective_basis_count=0
    nbasis=size(raw_basis,2)
    ok=size(raw_basis,1)>0 .and. nbasis>0 .and. &
      all(shape(orthonormal_basis)==shape(raw_basis)) .and. &
      all(shape(basis_transform)==[nbasis,nbasis]) .and. quadrature_weight>0d0 .and. &
      relative_cutoff>=100d0*epsilon(1d0) .and. relative_cutoff<1d0 .and. &
      ieee_is_finite(quadrature_weight) .and. ieee_is_finite(relative_cutoff) .and. &
      all(ieee_is_finite(real(raw_basis,8))) .and. all(ieee_is_finite(aimag(raw_basis)))
    if(.not.ok) then
      message='DG DC local basis: invalid fragment core basis'
      return
    end if
    maximum_basis_magnitude=maxval(abs(raw_basis))
    safe_basis_magnitude=sqrt(huge(1d0)/(quadrature_weight*real(size(raw_basis,1),8)))
    if(maximum_basis_magnitude>safe_basis_magnitude) then
      ok=.false.
      message='DG DC local basis: fragment core basis would overflow its metric'
      return
    end if
    lwork=max(1,2*nbasis-1)
    allocate(gram(nbasis,nbasis),work(lwork),eigenvalues(nbasis),rwork(max(1,3*nbasis-2)), &
      stat=allocation_status)
    if(allocation_status/=0) then
      ok=.false.; message='DG DC local basis: core metric workspace allocation failed'; return
    end if
    gram=quadrature_weight*matmul(conjg(transpose(raw_basis)),raw_basis)
    ok=all(ieee_is_finite(real(gram,8))) .and. all(ieee_is_finite(aimag(gram)))
    if(.not.ok) then
      message='DG DC local basis: nonfinite derived core metric'
      deallocate(gram,work,eigenvalues,rwork)
      return
    end if
    call zheev('V','U',nbasis,gram,nbasis,eigenvalues,work,lwork,rwork,info)
    ok=info==0
    if(ok) then
      maximum_eigenvalue=eigenvalues(nbasis)
      ok=ieee_is_finite(maximum_eigenvalue) .and. maximum_eigenvalue>0d0
    end if
    if(ok) then
      do source_column=nbasis,1,-1
        if(eigenvalues(source_column)<=relative_cutoff*maximum_eigenvalue) cycle
        effective_basis_count=effective_basis_count+1
        output_column=effective_basis_count
        basis_transform(:,output_column)=gram(:,source_column)/sqrt(eigenvalues(source_column))
      end do
      ok=effective_basis_count>0
    end if
    if(ok) then
      orthonormal_basis(:,1:effective_basis_count)= &
        matmul(raw_basis,basis_transform(:,1:effective_basis_count))
      ok=all(ieee_is_finite(real(orthonormal_basis,8))) .and. &
        all(ieee_is_finite(aimag(orthonormal_basis)))
      if(ok) then
        gram=(0d0,0d0)
        gram(1:effective_basis_count,1:effective_basis_count)=quadrature_weight* &
          matmul(conjg(transpose(orthonormal_basis(:,1:effective_basis_count))), &
          orthonormal_basis(:,1:effective_basis_count))
        do source_column=1,effective_basis_count
          gram(source_column,source_column)=gram(source_column,source_column)-1d0
        end do
        ok=maxval(abs(gram(1:effective_basis_count,1:effective_basis_count)))<=1d-10
      end if
    end if
    if(ok) then
      message=''
    else
      orthonormal_basis=(0d0,0d0)
      basis_transform=(0d0,0d0)
      effective_basis_count=0
      message='DG DC local basis: singular or invalid fragment core metric'
    end if
    deallocate(gram,work,eigenvalues,rwork)
  end subroutine orthonormalize_dg_dc_fragment_core_basis

  subroutine compose_dg_dc_distributed_hamiltonian_rows(layout,volume_block,interface_rows,lambda, &
      hamiltonian_rows,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(in) :: volume_block(:,:),interface_rows(:,:)
    real(8), intent(in) :: lambda
    complex(8), intent(out) :: hamiltonian_rows(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: rank,first,last

    hamiltonian_rows=(0d0,0d0)
    rank=0
#ifdef USE_MPI
    ! The layout arrays are rank ordered; fragment IDs may be permuted.
    do first=1,size(layout%fragment_ids)
      if(layout%fragment_ids(first)==layout%fragment_id) rank=first-1
    end do
#endif
    first=layout%basis_offsets(rank)+1
    last=layout%basis_offsets(rank+1)
    ok=layout%local_basis_count>0 .and. last-first+1==layout%local_basis_count .and. &
      all(shape(volume_block)==[layout%local_basis_count,layout%local_basis_count]) .and. &
      all(shape(interface_rows)==[layout%local_basis_count,layout%global_basis_count]) .and. &
      all(shape(hamiltonian_rows)==shape(interface_rows)) .and. lambda>=0d0 .and. lambda<=1d0 .and. &
      ieee_is_finite(lambda) .and. all(ieee_is_finite(real(volume_block,8))) .and. &
      all(ieee_is_finite(aimag(volume_block))) .and. &
      all(ieee_is_finite(real(interface_rows,8))) .and. all(ieee_is_finite(aimag(interface_rows)))
    if(.not.ok) then
      message='DG DC local basis: invalid distributed Hamiltonian rows'
      return
    end if
    hamiltonian_rows=lambda*interface_rows
    hamiltonian_rows(:,first:last)=hamiltonian_rows(:,first:last)+volume_block
    ok=all(ieee_is_finite(real(hamiltonian_rows,8))) .and. all(ieee_is_finite(aimag(hamiltonian_rows)))
    if(ok) then
      message=''
    else
      hamiltonian_rows=(0d0,0d0)
      message='DG DC local basis: nonfinite distributed Hamiltonian rows'
    end if
  end subroutine compose_dg_dc_distributed_hamiltonian_rows

  subroutine initialize_dg_dc_local_basis_layout(layout,fragment_id,local_basis_count,global_band_count, &
      geometry_fingerprint,operator_fingerprint,communicator,ok,message)
    type(s_dg_dc_local_basis_layout), intent(inout) :: layout
    integer, intent(in) :: fragment_id,local_basis_count,global_band_count,communicator
    integer(8), intent(in) :: geometry_fingerprint,operator_fingerprint
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: nproc,i,ierr,allocation_status
    integer(8) :: running_basis_count
    integer, allocatable :: basis_counts(:),band_counts(:)
    integer(8), allocatable :: geometry_values(:),operator_values(:)
    logical, allocatable :: seen_fragment(:)
    logical :: valid,stage_ok,all_stage_ok

    call release_dg_dc_local_basis_layout(layout)
    nproc=1
#ifdef USE_MPI
    call MPI_Comm_size(communicator,nproc,ierr)
    if(ierr/=MPI_SUCCESS) then
      ok=.false.; message='DG DC local basis: communicator size query failed'; return
    end if
#else
    if(communicator<0) then
      ok=.false.; message='DG DC local basis: invalid serial communicator'; return
    end if
#endif
    allocate(basis_counts(nproc),band_counts(nproc),layout%fragment_ids(nproc),geometry_values(nproc), &
      operator_values(nproc),layout%basis_offsets(0:nproc),seen_fragment(nproc),stat=allocation_status)
    stage_ok=allocation_status==0
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(.not.all_stage_ok) then
      if(allocated(basis_counts)) deallocate(basis_counts)
      if(allocated(band_counts)) deallocate(band_counts)
      if(allocated(geometry_values)) deallocate(geometry_values)
      if(allocated(operator_values)) deallocate(operator_values)
      if(allocated(seen_fragment)) deallocate(seen_fragment)
      call release_dg_dc_local_basis_layout(layout)
      ok=.false.; message='DG DC local basis: layout allocation failed collectively'; return
    end if
#ifdef USE_MPI
    call MPI_Allgather(local_basis_count,1,MPI_INTEGER,basis_counts,1,MPI_INTEGER,communicator,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(global_band_count,1,MPI_INTEGER,band_counts,1,MPI_INTEGER, &
      communicator,ierr)
    stage_ok=all_stage_ok .and. ierr==MPI_SUCCESS
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(fragment_id,1,MPI_INTEGER,layout%fragment_ids,1,MPI_INTEGER, &
      communicator,ierr)
    stage_ok=all_stage_ok .and. ierr==MPI_SUCCESS
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(geometry_fingerprint,1,MPI_INTEGER8,geometry_values,1, &
      MPI_INTEGER8,communicator,ierr)
    stage_ok=all_stage_ok .and. ierr==MPI_SUCCESS
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(all_stage_ok) call MPI_Allgather(operator_fingerprint,1,MPI_INTEGER8,operator_values,1, &
      MPI_INTEGER8,communicator,ierr)
    stage_ok=all_stage_ok .and. ierr==MPI_SUCCESS
    call collective_logical_and(stage_ok,communicator,all_stage_ok)
    if(.not.all_stage_ok) then
      deallocate(basis_counts,band_counts,geometry_values,operator_values,seen_fragment)
      call release_dg_dc_local_basis_layout(layout)
      ok=.false.; message='DG DC local basis: layout gather failed collectively'; return
    end if
#else
    basis_counts=[local_basis_count]
    band_counts=[global_band_count]
    layout%fragment_ids=[fragment_id]
    geometry_values=[geometry_fingerprint]
    operator_values=[operator_fingerprint]
#endif
    valid=all(basis_counts>0) .and. all(band_counts==band_counts(1)) .and. band_counts(1)>0 .and. &
      all(layout%fragment_ids>0) .and. all(layout%fragment_ids<=nproc) .and. &
      all(geometry_values==geometry_values(1)) .and. &
      all(operator_values==operator_values(1))
    seen_fragment=.false.
    do i=1,nproc
      if(layout%fragment_ids(i)>=1 .and. layout%fragment_ids(i)<=nproc) then
        valid=valid .and. .not.seen_fragment(layout%fragment_ids(i))
        seen_fragment(layout%fragment_ids(i))=.true.
      end if
    end do
    valid=valid .and. all(seen_fragment)
    layout%basis_offsets(0)=0
    running_basis_count=0_8
    do i=1,nproc
      if(int(basis_counts(i),8)>int(huge(1),8)-running_basis_count) then
        valid=.false.
        layout%basis_offsets(i)=0
      else
        running_basis_count=running_basis_count+int(basis_counts(i),8)
        layout%basis_offsets(i)=int(running_basis_count)
      end if
    end do
    if(.not.valid) then
      call release_dg_dc_local_basis_layout(layout)
      ok=.false.; message='DG DC local basis: invalid or rank-disagreeing layout'; return
    end if
    layout%fragment_id=fragment_id
    layout%local_basis_count=local_basis_count
    layout%global_basis_count=layout%basis_offsets(nproc)
    layout%global_band_count=global_band_count
    layout%geometry_fingerprint=geometry_fingerprint
    layout%operator_fingerprint=operator_fingerprint
    ok=.true.
    message=''
    deallocate(basis_counts,band_counts,geometry_values,operator_values,seen_fragment)
  end subroutine initialize_dg_dc_local_basis_layout

  subroutine release_dg_dc_local_basis_layout(layout)
    type(s_dg_dc_local_basis_layout), intent(inout) :: layout
    if(allocated(layout%basis_offsets)) deallocate(layout%basis_offsets)
    if(allocated(layout%fragment_ids)) deallocate(layout%fragment_ids)
    layout%fragment_id=0
    layout%local_basis_count=0
    layout%global_basis_count=0
    layout%global_band_count=0
    layout%geometry_fingerprint=0_8
    layout%operator_fingerprint=0_8
  end subroutine release_dg_dc_local_basis_layout

  subroutine assemble_dg_dc_local_basis_sipg_pair(local_value,local_outward_normal,neighbor_value, &
      neighbor_outward_normal,h_normal,face_weight,penalty_factor,physical_face,canonical_owner, &
      auxiliary_periodic_wrap,periodic_shift,expected_periodic_shift,neighbor_fragment,expected_neighbor_fragment, &
      h_local_local,h_local_neighbor,h_neighbor_local,h_neighbor_neighbor,ok,message)
    complex(8), intent(in) :: local_value(:,:),local_outward_normal(:,:)
    complex(8), intent(in) :: neighbor_value(:,:),neighbor_outward_normal(:,:)
    real(8), intent(in) :: h_normal,face_weight,penalty_factor
    logical, intent(in) :: physical_face,canonical_owner,auxiliary_periodic_wrap
    integer, intent(in) :: periodic_shift(3),expected_periodic_shift(3)
    integer, intent(in) :: neighbor_fragment,expected_neighbor_fragment
    complex(8), intent(out) :: h_local_local(:,:),h_local_neighbor(:,:)
    complex(8), intent(out) :: h_neighbor_local(:,:),h_neighbor_neighbor(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    type(s_dg_nodal_sipg_action) :: action
    integer :: ipoint,itrial,itest,info,nlocal,nneighbor
    complex(8) :: zero

    h_local_local=(0d0,0d0)
    h_local_neighbor=(0d0,0d0)
    h_neighbor_local=(0d0,0d0)
    h_neighbor_neighbor=(0d0,0d0)
    nlocal=size(local_value,2)
    nneighbor=size(neighbor_value,2)
    ok=size(local_value,1)>0 .and. nlocal>0 .and. nneighbor>0 .and. &
      all(shape(local_outward_normal)==shape(local_value)) .and. &
      size(neighbor_value,1)==size(local_value,1) .and. &
      all(shape(neighbor_outward_normal)==shape(neighbor_value)) .and. &
      all(shape(h_local_local)==[nlocal,nlocal]) .and. &
      all(shape(h_local_neighbor)==[nlocal,nneighbor]) .and. &
      all(shape(h_neighbor_local)==[nneighbor,nlocal]) .and. &
      all(shape(h_neighbor_neighbor)==[nneighbor,nneighbor]) .and. &
      h_normal>0d0 .and. face_weight>=0d0 .and. penalty_factor>0d0 .and. &
      ieee_is_finite(h_normal) .and. ieee_is_finite(face_weight) .and. ieee_is_finite(penalty_factor) .and. &
      all(abs(periodic_shift)<=1) .and. count(periodic_shift/=0)<=1 .and. &
      all(abs(expected_periodic_shift)<=1) .and. count(expected_periodic_shift/=0)<=1 .and. &
      neighbor_fragment>0 .and. expected_neighbor_fragment>0
    if(physical_face) ok=ok .and. neighbor_fragment==expected_neighbor_fragment .and. &
      all(periodic_shift==expected_periodic_shift)
    if(.not.ok) then
      message='DG DC local basis: invalid SIPG pair layout or controls'
      return
    end if
    if(.not.physical_face .or. .not.canonical_owner .or. auxiliary_periodic_wrap) then
      message=''
      return
    end if
    zero=(0d0,0d0)
    do ipoint=1,size(local_value,1)
      do itrial=1,nlocal
        call evaluate_dg_nodal_sipg_face(local_value(ipoint,itrial),zero, &
          local_outward_normal(ipoint,itrial),zero,h_normal,face_weight,penalty_factor,action,info)
        if(info/=0) ok=.false.
        do itest=1,nlocal
          h_local_local(itest,itrial)=h_local_local(itest,itrial)+ &
            conjg(local_value(ipoint,itest))*action%total_value(1)+ &
            conjg(local_outward_normal(ipoint,itest))*action%total_normal(1)
        end do
        do itest=1,nneighbor
          h_neighbor_local(itest,itrial)=h_neighbor_local(itest,itrial)+ &
            conjg(neighbor_value(ipoint,itest))*action%total_value(2)- &
            conjg(neighbor_outward_normal(ipoint,itest))*action%total_normal(2)
        end do
      end do
      do itrial=1,nneighbor
        call evaluate_dg_nodal_sipg_face(zero,neighbor_value(ipoint,itrial),zero, &
          -neighbor_outward_normal(ipoint,itrial),h_normal,face_weight,penalty_factor,action,info)
        if(info/=0) ok=.false.
        do itest=1,nlocal
          h_local_neighbor(itest,itrial)=h_local_neighbor(itest,itrial)+ &
            conjg(local_value(ipoint,itest))*action%total_value(1)+ &
            conjg(local_outward_normal(ipoint,itest))*action%total_normal(1)
        end do
        do itest=1,nneighbor
          h_neighbor_neighbor(itest,itrial)=h_neighbor_neighbor(itest,itrial)+ &
            conjg(neighbor_value(ipoint,itest))*action%total_value(2)- &
            conjg(neighbor_outward_normal(ipoint,itest))*action%total_normal(2)
        end do
      end do
    end do
    ok=ok .and. all(ieee_is_finite(real(h_local_local,8))) .and. &
      all(ieee_is_finite(aimag(h_local_local))) .and. &
      all(ieee_is_finite(real(h_local_neighbor,8))) .and. &
      all(ieee_is_finite(aimag(h_local_neighbor))) .and. &
      all(ieee_is_finite(real(h_neighbor_local,8))) .and. &
      all(ieee_is_finite(aimag(h_neighbor_local))) .and. &
      all(ieee_is_finite(real(h_neighbor_neighbor,8))) .and. &
      all(ieee_is_finite(aimag(h_neighbor_neighbor)))
    if(ok) then
      message=''
    else
      h_local_local=(0d0,0d0)
      h_local_neighbor=(0d0,0d0)
      h_neighbor_local=(0d0,0d0)
      h_neighbor_neighbor=(0d0,0d0)
      message='DG DC local basis: SIPG pair evaluation failed'
    end if
  end subroutine assemble_dg_dc_local_basis_sipg_pair

  subroutine compose_dg_dc_local_basis_pair_hamiltonian(volume_left,volume_right,h_local_local, &
      h_local_neighbor,h_neighbor_local,h_neighbor_neighbor,lambda,hamiltonian,ok,message)
    complex(8), intent(in) :: volume_left(:,:),volume_right(:,:),h_local_local(:,:),h_local_neighbor(:,:)
    complex(8), intent(in) :: h_neighbor_local(:,:),h_neighbor_neighbor(:,:)
    real(8), intent(in) :: lambda
    complex(8), intent(out) :: hamiltonian(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: nleft,nright
    nleft=size(volume_left,1)
    nright=size(volume_right,1)
    hamiltonian=(0d0,0d0)
    ok=nleft>0 .and. nright>0 .and. size(volume_left,2)==nleft .and. &
      size(volume_right,2)==nright .and. all(shape(h_local_local)==[nleft,nleft]) .and. &
      all(shape(h_local_neighbor)==[nleft,nright]) .and. &
      all(shape(h_neighbor_local)==[nright,nleft]) .and. &
      all(shape(h_neighbor_neighbor)==[nright,nright]) .and. &
      all(shape(hamiltonian)==[nleft+nright,nleft+nright]) .and. &
      ieee_is_finite(lambda) .and. lambda>=0d0 .and. lambda<=1d0
    if(.not.ok) then
      message='DG DC local basis: invalid pair Hamiltonian layout or lambda'
      return
    end if
    hamiltonian(1:nleft,1:nleft)=volume_left+lambda*h_local_local
    hamiltonian(1:nleft,nleft+1:nleft+nright)=lambda*h_local_neighbor
    hamiltonian(nleft+1:nleft+nright,1:nleft)=lambda*h_neighbor_local
    hamiltonian(nleft+1:nleft+nright,nleft+1:nleft+nright)=volume_right+lambda*h_neighbor_neighbor
    ok=all(ieee_is_finite(real(hamiltonian,8))) .and. all(ieee_is_finite(aimag(hamiltonian)))
    if(ok) then
      message=''
    else
      hamiltonian=(0d0,0d0)
      message='DG DC local basis: nonfinite pair Hamiltonian'
    end if
  end subroutine compose_dg_dc_local_basis_pair_hamiltonian

  subroutine solve_dg_dc_local_basis_bands_reference(layout,overlap,hamiltonian,communicator,metric_relative_tolerance, &
      eigenvalues,coefficients,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(in) :: overlap(:,:),hamiltonian(:,:)
    integer, intent(in) :: communicator
    real(8), intent(in) :: metric_relative_tolerance
    real(8), intent(out) :: eigenvalues(:)
    complex(8), intent(out) :: coefficients(:,:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: overlap_work(:,:),hamiltonian_work(:,:),reference_overlap(:,:), &
      reference_hamiltonian(:,:),work(:)
    real(8), allocatable :: all_eigenvalues(:),overlap_eigenvalues(:),rwork(:)
    integer :: nbasis,nband,lwork,info,allocation_status,ierr,rank
    real(8) :: hscale,sscale
    logical :: collective_ok,dimensions_agree,tolerance_agrees

    eigenvalues=0d0
    coefficients=(0d0,0d0)
    nbasis=layout%global_basis_count
    nband=layout%global_band_count
    call collective_integer_agreement(nbasis,communicator,dimensions_agree)
    call collective_integer_agreement(nband,communicator,collective_ok)
    dimensions_agree=dimensions_agree .and. collective_ok
    call collective_real_agreement(metric_relative_tolerance,communicator,tolerance_agrees)
    ok=nbasis>0 .and. nband>0 .and. nband<=nbasis .and. size(eigenvalues)==nband .and. &
      all(shape(overlap)==[nbasis,nbasis]) .and. all(shape(hamiltonian)==[nbasis,nbasis]) .and. &
      all(shape(coefficients)==[nbasis,nband]) .and. &
      metric_relative_tolerance>0d0 .and. metric_relative_tolerance<1d0 .and. &
      ieee_is_finite(metric_relative_tolerance) .and. &
      all(ieee_is_finite(real(overlap,8))) .and. all(ieee_is_finite(aimag(overlap))) .and. &
      all(ieee_is_finite(real(hamiltonian,8))) .and. all(ieee_is_finite(aimag(hamiltonian))) .and. &
      dimensions_agree .and. tolerance_agrees
    if(ok) then
      hscale=max(maxval(abs(hamiltonian)),tiny(1d0))
      sscale=max(maxval(abs(overlap)),tiny(1d0))
      ok=maxval(abs(hamiltonian-conjg(transpose(hamiltonian))))<=1d-12*hscale .and. &
        maxval(abs(overlap-conjg(transpose(overlap))))<=1d-12*sscale
    end if
    call collective_logical_and(ok,communicator,collective_ok)
    ok=collective_ok
    if(.not.ok) then
      message='DG DC local basis: invalid global eigensystem layout'
      return
    end if
    lwork=max(1,2*nbasis-1)
    allocate(overlap_work(nbasis,nbasis),hamiltonian_work(nbasis,nbasis), &
      reference_overlap(nbasis,nbasis),reference_hamiltonian(nbasis,nbasis),work(lwork), &
      all_eigenvalues(nbasis),overlap_eigenvalues(nbasis),rwork(max(1,3*nbasis-2)),stat=allocation_status)
    call collective_logical_and(allocation_status==0,communicator,collective_ok)
    if(.not.collective_ok) then
      ok=.false.
      message='DG DC local basis: eigensystem workspace allocation failed'
      return
    end if
    reference_overlap=overlap
    reference_hamiltonian=hamiltonian
    overlap_work=(0d0,0d0)
    hamiltonian_work=(0d0,0d0)
    all_eigenvalues=0d0
    overlap_eigenvalues=0d0
    rank=0
#ifdef USE_MPI
    call MPI_Comm_rank(communicator,rank,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,collective_ok)
    ok=collective_ok
    call MPI_Bcast(reference_overlap,nbasis*nbasis,MPI_DOUBLE_COMPLEX,0,communicator,ierr)
    call collective_logical_and(ok .and. ierr==MPI_SUCCESS,communicator,collective_ok)
    ok=collective_ok
    call MPI_Bcast(reference_hamiltonian,nbasis*nbasis,MPI_DOUBLE_COMPLEX,0,communicator,ierr)
    call collective_logical_and(ok .and. ierr==MPI_SUCCESS,communicator,collective_ok)
    ok=collective_ok
#endif
    ok=ok .and. maxval(abs(overlap-reference_overlap))<=1d-13*sscale .and. &
      maxval(abs(hamiltonian-reference_hamiltonian))<=1d-13*hscale
    call collective_logical_and(ok,communicator,collective_ok)
    ok=collective_ok
    if(ok .and. rank==0) then
      overlap_work=reference_overlap
      call zheev('N','U',nbasis,overlap_work,nbasis,overlap_eigenvalues,work,lwork,rwork,info)
      ok=info==0 .and. overlap_eigenvalues(1)>metric_relative_tolerance*overlap_eigenvalues(nbasis)
      if(ok) then
        overlap_work=reference_overlap
        hamiltonian_work=reference_hamiltonian
        call zhegv(1,'V','U',nbasis,hamiltonian_work,nbasis,overlap_work,nbasis,all_eigenvalues, &
          work,lwork,rwork,info)
        ok=info==0 .and. all(ieee_is_finite(all_eigenvalues)) .and. &
          all(ieee_is_finite(real(hamiltonian_work,8))) .and. all(ieee_is_finite(aimag(hamiltonian_work)))
      end if
    end if
#ifdef USE_MPI
    call MPI_Bcast(ok,1,MPI_LOGICAL,0,communicator,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,collective_ok)
    ok=ok .and. collective_ok
    call MPI_Bcast(all_eigenvalues,nbasis,MPI_DOUBLE_PRECISION,0,communicator,ierr)
    call collective_logical_and(ok .and. ierr==MPI_SUCCESS,communicator,collective_ok)
    ok=collective_ok
    call MPI_Bcast(hamiltonian_work,nbasis*nbasis,MPI_DOUBLE_COMPLEX,0,communicator,ierr)
    call collective_logical_and(ok .and. ierr==MPI_SUCCESS,communicator,collective_ok)
    ok=collective_ok
#endif
    if(ok) then
      eigenvalues=all_eigenvalues(1:nband)
      coefficients=hamiltonian_work(:,1:nband)
      message=''
    else
      message='DG DC local basis: generalized Hermitian eigensolve failed'
    end if
    deallocate(overlap_work,hamiltonian_work,reference_overlap,reference_hamiltonian,work, &
      all_eigenvalues,overlap_eigenvalues,rwork)
  end subroutine solve_dg_dc_local_basis_bands_reference

  subroutine solve_dg_dc_local_basis_bands_cg(layout,hamiltonian_rows,overlap_rows,communicator, &
      maximum_iterations,tolerance,coefficient_rows,eigenvalues,ok,message)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    complex(8), intent(in) :: hamiltonian_rows(:,:),overlap_rows(:,:)
    integer, intent(in) :: communicator,maximum_iterations
    real(8), intent(in) :: tolerance
    complex(8), intent(inout) :: coefficient_rows(:,:)
    real(8), intent(out) :: eigenvalues(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: global_vector(:),hvector(:),svector(:),residual(:),direction(:), &
      hdirection(:),sdirection(:)
    complex(8) :: reduced_h(2,2),reduced_s(2,2),work(8),overlap_value
    real(8) :: reduced_eigenvalues(2),rwork(4),residual_norm,denominator
    integer :: iband,jband,iteration,info,allocation_status
    logical :: stage_ok,operation_ok

    eigenvalues=0d0
    ok=layout%local_basis_count>0 .and. layout%global_basis_count>=layout%global_band_count .and. &
      maximum_iterations>0 .and. tolerance>0d0 .and. ieee_is_finite(tolerance) .and. &
      all(shape(hamiltonian_rows)==[layout%local_basis_count,layout%global_basis_count]) .and. &
      all(shape(overlap_rows)==shape(hamiltonian_rows)) .and. &
      all(shape(coefficient_rows)==[layout%local_basis_count,layout%global_band_count]) .and. &
      size(eigenvalues)==layout%global_band_count .and. &
      all(ieee_is_finite(real(hamiltonian_rows,8))) .and. all(ieee_is_finite(aimag(hamiltonian_rows))) .and. &
      all(ieee_is_finite(real(overlap_rows,8))) .and. all(ieee_is_finite(aimag(overlap_rows))) .and. &
      all(ieee_is_finite(real(coefficient_rows,8))) .and. all(ieee_is_finite(aimag(coefficient_rows)))
    call collective_logical_and(ok,communicator,stage_ok)
    if(stage_ok) call validate_coefficient_cg_metadata(layout,maximum_iterations,tolerance,hamiltonian_rows, &
      overlap_rows,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; message='DG DC coefficient CG: invalid distributed layout'; return
    end if
    allocate(global_vector(layout%global_basis_count),hvector(layout%local_basis_count), &
      svector(layout%local_basis_count),residual(layout%local_basis_count), &
      direction(layout%local_basis_count),hdirection(layout%local_basis_count), &
      sdirection(layout%local_basis_count),stat=allocation_status)
    call collective_logical_and(allocation_status==0,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; message='DG DC coefficient CG: workspace allocation failed'; return
    end if
    do iband=1,layout%global_band_count
      call orthonormalize_coefficient_vector(coefficient_rows(:,iband),coefficient_rows(:,1:iband-1), &
        overlap_rows,layout,communicator,stage_ok)
      if(.not.stage_ok) exit
      do iteration=1,maximum_iterations
        call gather_coefficient_vector(coefficient_rows(:,iband),layout,communicator,global_vector,stage_ok)
        if(.not.stage_ok) exit
        hvector=matmul(hamiltonian_rows,global_vector)
        svector=matmul(overlap_rows,global_vector)
        call distributed_inner_product(coefficient_rows(:,iband),hvector,communicator,reduced_h(1,1),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        if(.not.stage_ok) exit
        call distributed_inner_product(coefficient_rows(:,iband),svector,communicator,reduced_s(1,1),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        if(.not.stage_ok .or. real(reduced_s(1,1),8)<=0d0) then
          stage_ok=.false.; exit
        end if
        eigenvalues(iband)=real(reduced_h(1,1)/reduced_s(1,1),8)
        residual=hvector-eigenvalues(iband)*svector
        call distributed_inner_product(residual,residual,communicator,overlap_value,operation_ok)
        stage_ok=stage_ok .and. operation_ok
        residual_norm=sqrt(max(0d0,real(overlap_value,8)))
        if(.not.stage_ok .or. .not.ieee_is_finite(residual_norm)) then
          stage_ok=.false.; exit
        end if
        if(residual_norm<=tolerance) exit
        direction=-residual
        do jband=1,iband-1
          call gather_coefficient_vector(direction,layout,communicator,global_vector,stage_ok)
          if(.not.stage_ok) exit
          sdirection=matmul(overlap_rows,global_vector)
          call distributed_inner_product(coefficient_rows(:,jband),sdirection,communicator,overlap_value,operation_ok)
          stage_ok=stage_ok .and. operation_ok
          if(.not.stage_ok) exit
          direction=direction-coefficient_rows(:,jband)*overlap_value
        end do
        if(.not.stage_ok) exit
        call gather_coefficient_vector(direction,layout,communicator,global_vector,stage_ok)
        if(.not.stage_ok) exit
        hdirection=matmul(hamiltonian_rows,global_vector)
        sdirection=matmul(overlap_rows,global_vector)
        call distributed_inner_product(coefficient_rows(:,iband),hvector,communicator,reduced_h(1,1),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(coefficient_rows(:,iband),hdirection,communicator,reduced_h(1,2),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(direction,hvector,communicator,reduced_h(2,1),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(direction,hdirection,communicator,reduced_h(2,2),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(coefficient_rows(:,iband),svector,communicator,reduced_s(1,1),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(coefficient_rows(:,iband),sdirection,communicator,reduced_s(1,2),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(direction,svector,communicator,reduced_s(2,1),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        call distributed_inner_product(direction,sdirection,communicator,reduced_s(2,2),operation_ok)
        stage_ok=stage_ok .and. operation_ok
        if(.not.stage_ok) exit
        call zhegv(1,'V','U',2,reduced_h,2,reduced_s,2,reduced_eigenvalues,work,8,rwork,info)
        if(info/=0) then
          stage_ok=.false.; exit
        end if
        coefficient_rows(:,iband)=coefficient_rows(:,iband)*reduced_h(1,1)+direction*reduced_h(2,1)
        eigenvalues(iband)=reduced_eigenvalues(1)
      end do
      if(.not.stage_ok .or. iteration>maximum_iterations) then
        stage_ok=.false.; exit
      end if
    end do
    call collective_logical_and(stage_ok,communicator,ok)
    if(ok) then
      message=''
    else
      coefficient_rows=(0d0,0d0)
      eigenvalues=0d0
      message='DG DC coefficient CG: convergence or collective operation failed'
    end if
    deallocate(global_vector,hvector,svector,residual,direction,hdirection,sdirection)
  end subroutine solve_dg_dc_local_basis_bands_cg

  subroutine validate_coefficient_cg_metadata(layout,maximum_iterations,tolerance,hamiltonian_rows, &
      overlap_rows,communicator,ok)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    integer, intent(in) :: maximum_iterations,communicator
    real(8), intent(in) :: tolerance
    complex(8), intent(in) :: hamiltonian_rows(:,:),overlap_rows(:,:)
    logical, intent(out) :: ok
    integer :: nproc,rank,i,j,global_row,ierr
    integer, allocatable :: reference_offsets(:),reference_fragments(:)
    logical :: agrees,stage_ok
    real(8) :: expected
    complex(8), allocatable :: reference_row(:)
    real(8) :: local_defect,global_defect,hscale
#ifdef USE_MPI
    integer :: owner,local_row
#endif

    nproc=1; rank=0; ok=.true.
#ifdef USE_MPI
    call MPI_Comm_size(communicator,nproc,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,stage_ok)
    ok=ok .and. stage_ok
    call MPI_Comm_rank(communicator,rank,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,stage_ok)
    ok=ok .and. stage_ok
#endif
    call collective_integer_agreement(layout%global_basis_count,communicator,agrees); ok=ok .and. agrees
    call collective_integer_agreement(layout%global_band_count,communicator,agrees); ok=ok .and. agrees
    call collective_integer_agreement(maximum_iterations,communicator,agrees); ok=ok .and. agrees
    call collective_real_agreement(tolerance,communicator,agrees); ok=ok .and. agrees
    ok=ok .and. allocated(layout%basis_offsets) .and. allocated(layout%fragment_ids)
    if(ok) ok=size(layout%basis_offsets)==nproc+1 .and. size(layout%fragment_ids)==nproc
    call collective_logical_and(ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.; return
    end if
    allocate(reference_offsets(0:nproc),reference_fragments(nproc))
    reference_offsets=layout%basis_offsets
    reference_fragments=layout%fragment_ids
#ifdef USE_MPI
    call MPI_Bcast(reference_offsets,nproc+1,MPI_INTEGER,0,communicator,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,stage_ok); ok=ok .and. stage_ok
    call MPI_Bcast(reference_fragments,nproc,MPI_INTEGER,0,communicator,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,stage_ok); ok=ok .and. stage_ok
#endif
    ok=ok .and. all(layout%basis_offsets==reference_offsets) .and. &
      all(layout%fragment_ids==reference_fragments)
    ok=ok .and. reference_offsets(0)==0 .and. &
      reference_offsets(nproc)==layout%global_basis_count .and. &
      all(reference_offsets(1:nproc)>reference_offsets(0:nproc-1)) .and. &
      all(reference_offsets>=0) .and. all(reference_offsets<=layout%global_basis_count)
    if(ok) ok=layout%local_basis_count==reference_offsets(rank+1)-reference_offsets(rank)
    call collective_logical_and(ok,communicator,stage_ok)
    if(.not.stage_ok) then
      ok=.false.
      deallocate(reference_offsets,reference_fragments)
      return
    end if
    do j=1,layout%global_basis_count
      do i=1,layout%local_basis_count
        global_row=reference_offsets(rank)+i
        expected=merge(1d0,0d0,global_row==j)
        if(abs(overlap_rows(i,j)-expected)>1d-12) ok=.false.
      end do
    end do
    local_defect=0d0
    hscale=max(1d0,maxval(abs(hamiltonian_rows)))
#ifdef USE_MPI
    allocate(reference_row(layout%global_basis_count))
    do j=1,layout%global_basis_count
      owner=0
      do i=0,nproc-1
        if(j>reference_offsets(i) .and. j<=reference_offsets(i+1)) owner=i
      end do
      reference_row=(0d0,0d0)
      if(rank==owner) then
        local_row=j-reference_offsets(rank)
        reference_row=hamiltonian_rows(local_row,:)
      end if
      call MPI_Bcast(reference_row,layout%global_basis_count,MPI_DOUBLE_COMPLEX,owner,communicator,ierr)
      ok=ok .and. ierr==MPI_SUCCESS
      if(ierr==MPI_SUCCESS) then
        do i=1,layout%local_basis_count
          local_defect=max(local_defect,abs(hamiltonian_rows(i,j)- &
            conjg(reference_row(reference_offsets(rank)+i))))
        end do
      end if
    end do
    deallocate(reference_row)
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
    call MPI_Allreduce(hscale,expected,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
    hscale=expected
#else
    global_defect=maxval(abs(hamiltonian_rows-conjg(transpose(hamiltonian_rows))))
#endif
    ok=ok .and. global_defect<=1d-12*hscale
    call collective_logical_and(ok,communicator,stage_ok)
    ok=stage_ok
    deallocate(reference_offsets,reference_fragments)
  end subroutine validate_coefficient_cg_metadata

  subroutine orthonormalize_coefficient_vector(vector,previous,overlap_rows,layout,communicator,ok)
    complex(8), intent(inout) :: vector(:)
    complex(8), intent(in) :: previous(:,:)
    complex(8), intent(in) :: overlap_rows(:,:)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    complex(8), allocatable :: global_vector(:),svector(:)
    complex(8) :: value
    real(8) :: norm
    integer :: j
    allocate(global_vector(layout%global_basis_count),svector(layout%local_basis_count))
    ok=.true.
    do j=1,size(previous,2)
      call gather_coefficient_vector(vector,layout,communicator,global_vector,ok)
      if(.not.ok) return
      svector=matmul(overlap_rows,global_vector)
      call distributed_inner_product(previous(:,j),svector,communicator,value,ok)
      if(.not.ok) return
      vector=vector-previous(:,j)*value
    end do
    call gather_coefficient_vector(vector,layout,communicator,global_vector,ok)
    if(.not.ok) return
    svector=matmul(overlap_rows,global_vector)
    call distributed_inner_product(vector,svector,communicator,value,ok)
    norm=sqrt(max(0d0,real(value,8)))
    ok=ok .and. norm>1d-14 .and. ieee_is_finite(norm)
    if(ok) vector=vector/norm
    deallocate(global_vector,svector)
  end subroutine orthonormalize_coefficient_vector

  subroutine gather_coefficient_vector(local_vector,layout,communicator,global_vector,ok)
    complex(8), intent(in) :: local_vector(:)
    type(s_dg_dc_local_basis_layout), intent(in) :: layout
    integer, intent(in) :: communicator
    complex(8), intent(out) :: global_vector(:)
    logical, intent(out) :: ok
#ifdef USE_MPI
    integer :: ierr,nproc
    integer, allocatable :: counts(:),displacements(:)
    logical :: stage_ok
    call MPI_Comm_size(communicator,nproc,ierr)
    call collective_logical_and(ierr==MPI_SUCCESS,communicator,stage_ok)
    ok=stage_ok
    if(.not.ok) return
    ok=size(layout%basis_offsets)==nproc+1
    call collective_logical_and(ok,communicator,stage_ok)
    ok=stage_ok
    if(.not.ok) return
    allocate(counts(nproc),displacements(nproc))
    counts=layout%basis_offsets(1:nproc)-layout%basis_offsets(0:nproc-1)
    displacements=layout%basis_offsets(0:nproc-1)
    call MPI_Allgatherv(local_vector,size(local_vector),MPI_DOUBLE_COMPLEX,global_vector,counts, &
      displacements,MPI_DOUBLE_COMPLEX,communicator,ierr)
    ok=ierr==MPI_SUCCESS
    deallocate(counts,displacements)
#else
    ok=communicator>=0 .and. size(local_vector)==size(global_vector)
    if(ok) global_vector=local_vector
#endif
  end subroutine gather_coefficient_vector

  subroutine distributed_inner_product(left,right,communicator,value,ok)
    complex(8), intent(in) :: left(:),right(:)
    integer, intent(in) :: communicator
    complex(8), intent(out) :: value
    logical, intent(out) :: ok
    complex(8) :: local_value
#ifdef USE_MPI
    integer :: ierr
#endif
    ok=size(left)==size(right)
    local_value=(0d0,0d0)
    if(ok) local_value=sum(conjg(left)*right)
#ifdef USE_MPI
    call MPI_Allreduce(local_value,value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,communicator,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
#else
    value=local_value
    ok=ok .and. communicator>=0
#endif
  end subroutine distributed_inner_product

  subroutine reconstruct_dg_dc_local_basis_density(local_basis,local_coefficients,occupations,density,ok,message)
    complex(8), intent(in) :: local_basis(:,:),local_coefficients(:,:)
    real(8), intent(in) :: occupations(:)
    real(8), intent(out) :: density(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: local_states(:,:)
    integer :: iband,allocation_status
    density=0d0
    ok=size(local_basis,1)==size(density) .and. size(local_basis,2)>0 .and. &
      size(local_coefficients,1)==size(local_basis,2) .and. &
      size(local_coefficients,2)==size(occupations) .and. size(occupations)>0 .and. &
      all(occupations>=0d0) .and. all(ieee_is_finite(occupations)) .and. &
      all(ieee_is_finite(real(local_basis,8))) .and. all(ieee_is_finite(aimag(local_basis))) .and. &
      all(ieee_is_finite(real(local_coefficients,8))) .and. all(ieee_is_finite(aimag(local_coefficients)))
    if(.not.ok) then
      message='DG DC local basis: invalid density reconstruction layout'
      return
    end if
    allocate(local_states(size(local_basis,1),size(occupations)),stat=allocation_status)
    if(allocation_status/=0) then
      ok=.false.; message='DG DC local basis: density workspace allocation failed'; return
    end if
    local_states=matmul(local_basis,local_coefficients)
    do iband=1,size(occupations)
      density=density+occupations(iband)*abs(local_states(:,iband))**2
    end do
    ok=all(ieee_is_finite(density))
    if(ok) then
      message=''
    else
      density=0d0
      message='DG DC local basis: nonfinite reconstructed density'
    end if
    deallocate(local_states)
  end subroutine reconstruct_dg_dc_local_basis_density

  subroutine validate_dg_dc_local_basis_density(occupations,maximum_occupation,expected_electrons, &
      local_density,quadrature_weights,communicator,ok,message)
    real(8), intent(in) :: occupations(:),maximum_occupation,expected_electrons
    real(8), intent(in) :: local_density(:),quadrature_weights(:)
    integer, intent(in) :: communicator
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), allocatable :: reference_occupations(:)
    real(8) :: local_charge,global_charge,occupation_scale
    integer :: ierr,allocation_status
    logical :: collective_ok,occupation_size_agrees,maximum_agrees,electron_count_agrees

    call collective_integer_agreement(size(occupations),communicator,occupation_size_agrees)
    call collective_real_agreement(maximum_occupation,communicator,maximum_agrees)
    call collective_real_agreement(expected_electrons,communicator,electron_count_agrees)
    allocate(reference_occupations(size(occupations)),stat=allocation_status)
    ok=allocation_status==0 .and. size(occupations)>0 .and. size(local_density)>0 .and. &
      size(quadrature_weights)==size(local_density) .and. maximum_occupation>0d0 .and. &
      expected_electrons>=0d0 .and. all(occupations>=0d0) .and. &
      all(occupations<=maximum_occupation) .and. all(local_density>=0d0) .and. &
      all(quadrature_weights>=0d0) .and. all(ieee_is_finite(occupations)) .and. &
      all(ieee_is_finite(local_density)) .and. all(ieee_is_finite(quadrature_weights)) .and. &
      occupation_size_agrees .and. maximum_agrees .and. electron_count_agrees
    call collective_logical_and(ok,communicator,collective_ok)
    if(.not.collective_ok) then
      if(allocated(reference_occupations)) deallocate(reference_occupations)
      ok=.false.; message='DG DC local basis: invalid density normalization inputs'; return
    end if
    reference_occupations=occupations
#ifdef USE_MPI
    call MPI_Bcast(reference_occupations,size(occupations),MPI_DOUBLE_PRECISION,0,communicator,ierr)
    ok=ierr==MPI_SUCCESS
#else
    ok=communicator>=0
#endif
    occupation_scale=max(maximum_occupation,1d0)
    ok=ok .and. maxval(abs(occupations-reference_occupations))<=1d-13*occupation_scale .and. &
      abs(sum(occupations)-expected_electrons)<=1d-12*max(expected_electrons,1d0)
    local_charge=sum(local_density*quadrature_weights)
    global_charge=local_charge
#ifdef USE_MPI
    call MPI_Allreduce(local_charge,global_charge,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)
    ok=ok .and. ierr==MPI_SUCCESS
#endif
    ok=ok .and. abs(global_charge-expected_electrons)<=1d-10*max(expected_electrons,1d0)
    call collective_logical_and(ok,communicator,collective_ok)
    ok=collective_ok
    if(ok) then
      message=''
    else
      message='DG DC local basis: occupation or integrated charge mismatch'
    end if
    deallocate(reference_occupations)
  end subroutine validate_dg_dc_local_basis_density

  subroutine collective_integer_agreement(local_value,communicator,agrees)
    integer, intent(in) :: local_value,communicator
    logical, intent(out) :: agrees
    integer :: minimum_value,maximum_value,ierr
#ifdef USE_MPI
    call MPI_Allreduce(local_value,minimum_value,1,MPI_INTEGER,MPI_MIN,communicator,ierr)
    agrees=ierr==MPI_SUCCESS
    call MPI_Allreduce(local_value,maximum_value,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    agrees=agrees .and. ierr==MPI_SUCCESS .and. minimum_value==maximum_value
#else
    agrees=communicator>=0
#endif
  end subroutine collective_integer_agreement

  subroutine collective_real_agreement(local_value,communicator,agrees)
    real(8), intent(in) :: local_value
    integer, intent(in) :: communicator
    logical, intent(out) :: agrees
    real(8) :: minimum_value,maximum_value
    integer :: ierr
#ifdef USE_MPI
    call MPI_Allreduce(local_value,minimum_value,1,MPI_DOUBLE_PRECISION,MPI_MIN,communicator,ierr)
    agrees=ierr==MPI_SUCCESS
    call MPI_Allreduce(local_value,maximum_value,1,MPI_DOUBLE_PRECISION,MPI_MAX,communicator,ierr)
    agrees=agrees .and. ierr==MPI_SUCCESS .and. minimum_value==maximum_value
#else
    agrees=communicator>=0
#endif
  end subroutine collective_real_agreement

  subroutine collective_logical_and(local_value,communicator,global_value)
    logical, intent(in) :: local_value
    integer, intent(in) :: communicator
    logical, intent(out) :: global_value
#ifdef USE_MPI
    integer :: ierr
    call MPI_Allreduce(local_value,global_value,1,MPI_LOGICAL,MPI_LAND,communicator,ierr)
    if(ierr/=MPI_SUCCESS) global_value=.false.
#else
    global_value=local_value
    if(communicator<0) global_value=.false.
#endif
  end subroutine collective_logical_and
end module dg_dc_local_basis_ground_state
