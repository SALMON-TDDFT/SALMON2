!
! Matrix-free SIPG face action used directly by the DC orbital CG.
!
#include "config.h"
module dg_dc_direct_sipg
#ifdef USE_MPI
  use mpi
#endif
  use dg_nodal_sipg, only: s_dg_nodal_sipg_action,evaluate_dg_nodal_sipg_face
  use dg_buffer_window_projector, only: s_dg_buffer_projector_diagnostics, &
    build_dg_buffer_window_projector_dual
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: apply_dg_dc_direct_face_mpi,evaluate_dg_dc_direct_local_face
  public :: apply_dg_dc_frozen_projected_face
  public :: s_dg_dc_frozen_face_projector,freeze_dg_dc_face_projector
  public :: reconstruct_dg_dc_frozen_face
  public :: dg_dc_canonical_to_raw_index
  public :: apply_dg_dc_frozen_projected_face_plane

  type :: s_dg_dc_frozen_face_projector
    integer :: generation=0
    integer(8) :: operator_fingerprint=0_8
    real(8),allocatable :: core_basis(:,:),normal_basis(:,:),dual(:,:)
    type(s_dg_buffer_projector_diagnostics) :: diagnostics
  end type
contains
  pure integer function dg_dc_canonical_to_raw_index(index,core_count,buffer_count)
    integer,intent(in)::index,core_count,buffer_count
    if(index<=buffer_count)then
      dg_dc_canonical_to_raw_index=core_count+buffer_count+index
    else
      dg_dc_canonical_to_raw_index=index-buffer_count
    endif
  end function dg_dc_canonical_to_raw_index

  subroutine freeze_dg_dc_face_projector(neighbor_buffer_basis,neighbor_core_basis,&
      neighbor_core_normal,weights,rank_tolerance,generation,operator_fingerprint,&
      projector,ok,message)
    real(8),intent(in) :: neighbor_buffer_basis(:,:),neighbor_core_basis(:,:)
    real(8),intent(in) :: neighbor_core_normal(:,:),weights(:),rank_tolerance
    integer,intent(in) :: generation
    integer(8),intent(in) :: operator_fingerprint
    type(s_dg_dc_frozen_face_projector),intent(out) :: projector
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: npoint,nstate,allocation_status

    ok=.false.;message=''
    npoint=size(neighbor_buffer_basis,1);nstate=size(neighbor_buffer_basis,2)
    if(npoint<1.or.nstate<1.or.size(weights)/=npoint.or.&
       size(neighbor_core_basis,2)/=nstate.or.size(neighbor_core_normal,2)/=nstate.or.&
       size(neighbor_core_basis,1)/=size(neighbor_core_normal,1).or.&
       generation<=0.or.operator_fingerprint==0_8.or.rank_tolerance<0d0)then
      message='direct DC SIPG: invalid frozen face projector context';return
    endif
    allocate(projector%core_basis(size(neighbor_core_basis,1),nstate),&
      projector%normal_basis(size(neighbor_core_normal,1),nstate),&
      projector%dual(nstate,npoint),stat=allocation_status)
    if(allocation_status/=0)then
      ok=.false.;message='direct DC SIPG: frozen face operator allocation failed';return
    endif
    call build_dg_buffer_window_projector_dual(neighbor_buffer_basis,weights,rank_tolerance,&
      projector%dual,projector%diagnostics,ok,message)
    if(.not.ok)return
    projector%core_basis=neighbor_core_basis;projector%normal_basis=neighbor_core_normal
    projector%generation=generation
    projector%operator_fingerprint=operator_fingerprint
    message=''
  end subroutine freeze_dg_dc_face_projector

  subroutine reconstruct_dg_dc_frozen_face(projector,buffer_values,expected_generation,&
      expected_fingerprint,reconstructed_values,reconstructed_normals,ok,message)
    type(s_dg_dc_frozen_face_projector),intent(in) :: projector
    real(8),intent(in) :: buffer_values(:,:)
    integer,intent(in) :: expected_generation
    integer(8),intent(in) :: expected_fingerprint
    real(8),intent(out) :: reconstructed_values(:,:),reconstructed_normals(:,:)
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    real(8),allocatable :: coefficients(:,:)
    integer :: allocation_status
    ok=allocated(projector%core_basis).and.allocated(projector%normal_basis).and.&
      allocated(projector%dual).and.&
      projector%generation>0.and.projector%generation==expected_generation.and.&
      projector%operator_fingerprint/=0_8.and.projector%operator_fingerprint==expected_fingerprint
    if(ok)ok=size(buffer_values,1)==size(projector%dual,2).and.&
      size(reconstructed_values,1)==size(projector%core_basis,1).and.&
      size(reconstructed_values,2)==size(buffer_values,2).and.&
      all(shape(reconstructed_values)==shape(reconstructed_normals))
    reconstructed_values=0d0;reconstructed_normals=0d0
    if(.not.ok)then
      message='direct DC SIPG: stale generation or invalid frozen reconstruction';return
    endif
    allocate(coefficients(size(projector%dual,1),size(buffer_values,2)),stat=allocation_status)
    if(allocation_status/=0)then
      ok=.false.;message='direct DC SIPG: frozen reconstruction allocation failed';return
    endif
    coefficients=matmul(projector%dual,buffer_values)
    reconstructed_values=matmul(projector%core_basis,coefficients)
    reconstructed_normals=matmul(projector%normal_basis,coefficients)
    ok=all(ieee_is_finite(reconstructed_values)).and.all(ieee_is_finite(reconstructed_normals))
    if(ok)then
      message=''
    else
      reconstructed_values=0d0;reconstructed_normals=0d0
      message='direct DC SIPG: nonfinite frozen reconstruction'
    endif
  end subroutine reconstruct_dg_dc_frozen_face

  subroutine apply_dg_dc_frozen_projected_face_plane(projector,buffer_values,core_lines,&
      derivative_weights,expected_generation,expected_fingerprint,side,h_normal,face_weight,&
      penalty_factor,value_action,normal_lift,ok,message)
    type(s_dg_dc_frozen_face_projector),intent(in)::projector
    real(8),intent(in)::buffer_values(:,:),core_lines(:,:,:),derivative_weights(:)
    integer,intent(in)::expected_generation,side
    integer(8),intent(in)::expected_fingerprint
    real(8),intent(in)::h_normal,face_weight,penalty_factor
    real(8),intent(out)::value_action(:,:),normal_lift(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::neighbor_values(:,:),neighbor_normals(:,:)
    complex(8)::face_value,face_normal,local_normal
    integer::point,state,index,info,allocation_status
    value_action=0d0;normal_lift=0d0;ok=.false.;message=''
    if(size(core_lines,2)/=size(derivative_weights).or.&
      size(core_lines,3)/=size(buffer_values,2).or.&
      any(shape(value_action)/=[size(core_lines,1),size(core_lines,3)]).or.&
      any(shape(normal_lift)/=shape(core_lines)))then
      message='direct DC SIPG: invalid frozen face-plane dimensions';return
    endif
    allocate(neighbor_values(size(core_lines,1),size(core_lines,3)),&
      neighbor_normals(size(core_lines,1),size(core_lines,3)),stat=allocation_status)
    if(allocation_status/=0)then;message='direct DC SIPG: face-plane allocation failed';return;endif
    call reconstruct_dg_dc_frozen_face(projector,buffer_values,expected_generation,&
      expected_fingerprint,neighbor_values,neighbor_normals,ok,message)
    if(.not.ok)return
    do state=1,size(core_lines,3);do point=1,size(core_lines,1)
      local_normal=(0d0,0d0)
      do index=1,size(derivative_weights)
        local_normal=local_normal+derivative_weights(index)*core_lines(point,index,state)
      enddo
      call evaluate_dg_dc_direct_local_face(cmplx(core_lines(point,1,state),0d0,8),local_normal,&
        cmplx(neighbor_values(point,state),0d0,8),cmplx(neighbor_normals(point,state),0d0,8),&
        side,h_normal,face_weight,penalty_factor,1d0,face_value,face_normal,info)
      if(info/=0)then;ok=.false.;message='direct DC SIPG: face-plane evaluation failed';return;endif
      value_action(point,state)=real(face_value,8)
      do index=1,size(derivative_weights)
        normal_lift(point,index,state)=derivative_weights(index)*real(face_normal,8)
      enddo
    enddo;enddo
    ok=.true.;message=''
  end subroutine apply_dg_dc_frozen_projected_face_plane

  subroutine apply_dg_dc_frozen_projected_face(local_value,local_outward_normal,&
      reconstructed_neighbor_value,reconstructed_neighbor_normal,projector_generation,&
      expected_generation,side,h_normal,face_weight,penalty_factor,scale,&
      local_face_value,local_face_normal,ok,message)
    complex(8),intent(in) :: local_value(:,:,:,:),local_outward_normal(:,:,:,:)
    complex(8),intent(in) :: reconstructed_neighbor_value(:,:,:,:),reconstructed_neighbor_normal(:,:,:,:)
    integer,intent(in) :: projector_generation,expected_generation,side
    real(8),intent(in) :: h_normal,face_weight,penalty_factor,scale
    complex(8),intent(out) :: local_face_value(:,:,:,:),local_face_normal(:,:,:,:)
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    integer :: i,j,istate,ispin,info
    ok=projector_generation>0.and.projector_generation==expected_generation.and.abs(side)==1.and.&
      h_normal>0d0.and.face_weight>=0d0.and.penalty_factor>0d0.and.scale>=0d0.and.&
      all(shape(local_value)==shape(local_outward_normal)).and.&
      all(shape(local_value)==shape(reconstructed_neighbor_value)).and.&
      all(shape(local_value)==shape(reconstructed_neighbor_normal)).and.&
      all(shape(local_value)==shape(local_face_value)).and.all(shape(local_value)==shape(local_face_normal))
    local_face_value=(0d0,0d0);local_face_normal=(0d0,0d0)
    if(.not.ok)then
      message='direct DC SIPG: stale generation or invalid frozen projected face';return
    endif
    do ispin=1,size(local_value,4);do istate=1,size(local_value,3)
      do j=1,size(local_value,2);do i=1,size(local_value,1)
        call evaluate_dg_dc_direct_local_face(local_value(i,j,istate,ispin),&
          local_outward_normal(i,j,istate,ispin),reconstructed_neighbor_value(i,j,istate,ispin),&
          reconstructed_neighbor_normal(i,j,istate,ispin),side,h_normal,face_weight,penalty_factor,scale,&
          local_face_value(i,j,istate,ispin),local_face_normal(i,j,istate,ispin),info)
        if(info/=0)ok=.false.
      enddo;enddo
    enddo;enddo
    if(ok)then
      message=''
    else
      local_face_value=(0d0,0d0);local_face_normal=(0d0,0d0)
      message='direct DC SIPG: frozen projected face evaluation failed'
    endif
  end subroutine apply_dg_dc_frozen_projected_face

  pure subroutine evaluate_dg_dc_direct_local_face(local_value,local_outward_normal,neighbor_value, &
      neighbor_outward_normal,side,h_normal,face_weight,penalty_factor,scale,face_value,face_normal,info)
    complex(8),intent(in) :: local_value,local_outward_normal,neighbor_value,neighbor_outward_normal
    integer,intent(in) :: side
    real(8),intent(in) :: h_normal,face_weight,penalty_factor,scale
    complex(8),intent(out) :: face_value,face_normal
    integer,intent(out) :: info
    type(s_dg_nodal_sipg_action) :: action
    info=1
    face_value=(0d0,0d0)
    face_normal=(0d0,0d0)
    if(abs(side)/=1 .or. scale<0d0 .or. .not.ieee_is_finite(scale))return
    if(scale<=0d0)then
      info=0
      return
    end if
    if(side>0)then
      call evaluate_dg_nodal_sipg_face(local_value,neighbor_value,local_outward_normal, &
        -neighbor_outward_normal,h_normal,face_weight,penalty_factor,action,info)
      if(info==0)then
        face_value=scale*action%total_value(1)
        face_normal=scale*action%total_normal(1)
      end if
    else
      call evaluate_dg_nodal_sipg_face(neighbor_value,local_value,neighbor_outward_normal, &
        -local_outward_normal,h_normal,face_weight,penalty_factor,action,info)
      if(info==0)then
        face_value=scale*action%total_value(2)
        face_normal=-scale*action%total_normal(2)
      end if
    end if
  end subroutine evaluate_dg_dc_direct_local_face

  subroutine apply_dg_dc_direct_face_mpi(local_value,local_outward_normal,neighbor_rank,side, &
      h_normal,face_weight,penalty_factor,scale,communicator,local_face_value,local_face_normal,ok,message)
    complex(8),intent(in) :: local_value(:,:,:,:),local_outward_normal(:,:,:,:)
    integer,intent(in) :: neighbor_rank,side,communicator
    real(8),intent(in) :: h_normal,face_weight,penalty_factor,scale
    complex(8),intent(out) :: local_face_value(:,:,:,:),local_face_normal(:,:,:,:)
    logical,intent(out) :: ok
    character(*),intent(out) :: message
    complex(8),allocatable :: neighbor_value(:,:,:,:),neighbor_normal(:,:,:,:)
    integer :: i,j,istate,ispin,info
#ifdef USE_MPI
    integer :: ierr_value,ierr_normal,ierr_agree,local_valid,global_valid
#endif
    ok=abs(side)==1 .and. neighbor_rank>=0 .and. h_normal>0d0 .and. face_weight>=0d0 .and. &
      penalty_factor>0d0 .and. scale>=0d0 .and. ieee_is_finite(scale) .and. &
      all(shape(local_value)==shape(local_outward_normal)) .and. &
      all(shape(local_value)==shape(local_face_value)) .and. &
      all(shape(local_value)==shape(local_face_normal))
#ifdef USE_MPI
    local_valid=merge(1,0,ok)
    call MPI_Allreduce(local_valid,global_valid,1,MPI_INTEGER,MPI_MIN,communicator,ierr_agree)
    ok=ierr_agree==MPI_SUCCESS .and. global_valid==1
#endif
    if(.not.ok) then
      local_face_value=(0d0,0d0); local_face_normal=(0d0,0d0)
      message='direct DC SIPG: invalid face context collectively'
      return
    end if
    allocate(neighbor_value,mold=local_value)
    allocate(neighbor_normal,mold=local_outward_normal)
#ifdef USE_MPI
    call MPI_Sendrecv(local_value,size(local_value),MPI_DOUBLE_COMPLEX,neighbor_rank,971, &
      neighbor_value,size(neighbor_value),MPI_DOUBLE_COMPLEX,neighbor_rank,971,communicator, &
      MPI_STATUS_IGNORE,ierr_value)
    call MPI_Sendrecv(local_outward_normal,size(local_outward_normal),MPI_DOUBLE_COMPLEX,neighbor_rank,972, &
      neighbor_normal,size(neighbor_normal),MPI_DOUBLE_COMPLEX,neighbor_rank,972,communicator, &
      MPI_STATUS_IGNORE,ierr_normal)
    local_valid=merge(1,0,ierr_value==MPI_SUCCESS .and. ierr_normal==MPI_SUCCESS)
    call MPI_Allreduce(local_valid,global_valid,1,MPI_INTEGER,MPI_MIN,communicator,ierr_agree)
    ok=ierr_agree==MPI_SUCCESS .and. global_valid==1
#else
    ok=neighbor_rank==0 .and. communicator>=0
    neighbor_value=local_value
    neighbor_normal=local_outward_normal
#endif
    local_face_value=(0d0,0d0)
    local_face_normal=(0d0,0d0)
    if(ok .and. scale>0d0) then
      do ispin=1,size(local_value,4)
        do istate=1,size(local_value,3)
          do j=1,size(local_value,2)
            do i=1,size(local_value,1)
              call evaluate_dg_dc_direct_local_face(local_value(i,j,istate,ispin), &
                local_outward_normal(i,j,istate,ispin),neighbor_value(i,j,istate,ispin), &
                neighbor_normal(i,j,istate,ispin),side,h_normal,face_weight,penalty_factor,scale, &
                local_face_value(i,j,istate,ispin),local_face_normal(i,j,istate,ispin),info)
              ok=ok .and. info==0
            end do
          end do
        end do
      end do
    end if
    deallocate(neighbor_value,neighbor_normal)
    if(ok) then
      message=''
    else
      local_face_value=(0d0,0d0); local_face_normal=(0d0,0d0)
      message='direct DC SIPG: face exchange or evaluation failed collectively'
    end if
  end subroutine apply_dg_dc_direct_face_mpi
end module dg_dc_direct_sipg
