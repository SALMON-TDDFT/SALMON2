#include "config.h"
program test_dg_dc_direct_cg_mpi
#ifdef USE_MPI
  use mpi
#endif
  use dg_dc_direct_sipg, only: apply_dg_dc_direct_face_mpi,apply_dg_dc_frozen_projected_face,&
    s_dg_dc_frozen_face_projector,freeze_dg_dc_face_projector,reconstruct_dg_dc_frozen_face,&
    apply_dg_dc_frozen_projected_face_plane
  use dg_dc_direct_sipg, only: dg_dc_canonical_to_raw_index
  implicit none
  integer :: communicator,rank,nproc,ierr
  integer :: index
  integer :: expected_raw(8)
  complex(8) :: value(2,1,2,1),normal(2,1,2,1)
  complex(8) :: action_value_0(2,1,2,1),action_normal_0(2,1,2,1)
  complex(8) :: action_value_half(2,1,2,1),action_normal_half(2,1,2,1)
  complex(8) :: action_value_one(2,1,2,1),action_normal_one(2,1,2,1)
  complex(8) :: raw_neighbor_value(2,1,2,1),raw_neighbor_normal(2,1,2,1)
  complex(8) :: projected_neighbor_value(2,1,2,1),projected_neighbor_normal(2,1,2,1)
  complex(8) :: projected_action_value(2,1,2,1),projected_action_normal(2,1,2,1)
  complex(8) :: raw_action_value(2,1,2,1),raw_action_normal(2,1,2,1)
  logical :: ok
  character(256) :: message
  type(s_dg_dc_frozen_face_projector) :: frozen,rotated_frozen
  real(8) :: basis(3,2),rotated_basis(3,2),core_normal(3,2),rotated_normal(3,2),rotation(2,2)
  real(8) :: weights(3),vectors(3,3),reconstructed(3,3),reconstructed_normal(3,3)
  real(8) :: rotated_reconstructed(3,3),rotated_reconstructed_normal(3,3)
  real(8) :: core_lines(3,2,3),derivative_weights(2)
  real(8) :: plane_value(3,3),plane_lift(3,2,3),rotated_plane_value(3,3),rotated_plane_lift(3,2,3)
#ifdef USE_MPI
  call MPI_Init(ierr)
  communicator=MPI_COMM_WORLD
  call MPI_Comm_rank(communicator,rank,ierr)
  call MPI_Comm_size(communicator,nproc,ierr)
#else
  communicator=0; rank=0; nproc=1
#endif
  if(nproc/=2) error stop 'fixture requires two ranks'
  expected_raw=[7,8,1,2,3,4,5,6]
  do index=1,8
    if(dg_dc_canonical_to_raw_index(index,4,2)/=expected_raw(index))&
      error stop 'canonical-to-raw DC storage mapping failed'
  enddo
  basis=reshape([1d0,0.25d0,-0.5d0,-0.5d0,1.5d0,0.75d0],shape(basis))
  core_normal=reshape([0.2d0,-0.1d0,0.3d0,0.7d0,0.4d0,-0.2d0],shape(core_normal))
  rotation=reshape([sqrt(0.5d0),sqrt(0.5d0),-sqrt(0.5d0),sqrt(0.5d0)],shape(rotation))
  rotated_basis=matmul(basis,rotation)
  rotated_normal=matmul(core_normal,rotation)
  weights=[0.2d0,0.3d0,0.5d0]
  vectors=reshape([1d0,2d0,-3d0,0.5d0,4d0,-2d0,-1d0,0.25d0,3d0],shape(vectors))
  call freeze_dg_dc_face_projector(basis,basis,core_normal,weights,1d-12,17,901_8,frozen,ok,message)
  if(.not.ok)error stop trim(message)
  call freeze_dg_dc_face_projector(rotated_basis,rotated_basis,rotated_normal,weights,1d-12,17,901_8,&
    rotated_frozen,ok,message)
  if(.not.ok)error stop trim(message)
  call reconstruct_dg_dc_frozen_face(frozen,vectors,17,901_8,reconstructed,reconstructed_normal,ok,message)
  if(.not.ok)error stop trim(message)
  call reconstruct_dg_dc_frozen_face(rotated_frozen,vectors,17,901_8,rotated_reconstructed,&
    rotated_reconstructed_normal,ok,message)
  if(.not.ok)error stop trim(message)
  if(maxval(abs(reconstructed-rotated_reconstructed))>1d-12.or.&
     maxval(abs(reconstructed_normal-rotated_reconstructed_normal))>1d-12)&
    error stop 'frozen face operator is not invariant to state-window rotation'
  core_lines(:,1,:)=reshape([1d0,2d0,3d0,0.5d0,-1d0,4d0,2d0,0.25d0,-2d0],[3,3])
  core_lines(:,2,:)=0.75d0*core_lines(:,1,:)
  derivative_weights=[-1.5d0,1.5d0]
  call apply_dg_dc_frozen_projected_face_plane(frozen,vectors,core_lines,derivative_weights,&
    17,901_8,1,0.4d0,0.16d0,9d0,plane_value,plane_lift,ok,message)
  if(.not.ok)error stop trim(message)
  call apply_dg_dc_frozen_projected_face_plane(rotated_frozen,vectors,core_lines,derivative_weights,&
    17,901_8,1,0.4d0,0.16d0,9d0,rotated_plane_value,rotated_plane_lift,ok,message)
  if(.not.ok)error stop trim(message)
  if(maxval(abs(plane_value-rotated_plane_value))>1d-12.or.&
     maxval(abs(plane_lift-rotated_plane_lift))>1d-12)&
    error stop 'production face-plane kernel is not invariant to state-window rotation'
  call reconstruct_dg_dc_frozen_face(frozen,vectors,18,901_8,reconstructed,reconstructed_normal,ok,message)
  if(ok)error stop 'stale reusable face projector generation was accepted'
  value=cmplx(dble(rank+1),0d0,8)
  value(:,:,2,:)=2d0*value(:,:,1,:)
  normal=cmplx(0.25d0*dble(2*rank-1),0d0,8)
  normal(:,:,2,:)=3d0*normal(:,:,1,:)
  call apply_dg_dc_direct_face_mpi(value,normal,1-rank,merge(1,-1,rank==0),0.4d0,0.16d0,9d0,0d0, &
    communicator,action_value_0,action_normal_0,ok,message)
  if(.not.ok) error stop trim(message)
  call apply_dg_dc_direct_face_mpi(value,normal,1-rank,merge(1,-1,rank==0),0.4d0,0.16d0,9d0,0.5d0, &
    communicator,action_value_half,action_normal_half,ok,message)
  if(.not.ok) error stop trim(message)
  call apply_dg_dc_direct_face_mpi(value,normal,1-rank,merge(1,-1,rank==0),0.4d0,0.16d0,9d0,1d0, &
    communicator,action_value_one,action_normal_one,ok,message)
  if(.not.ok) error stop trim(message)
  if(maxval(abs(action_value_0))+maxval(abs(action_normal_0))>1d-14) error stop 'scale zero is not zero'
  if(maxval(abs(2d0*action_value_half-action_value_one))>1d-13) error stop 'value action is not linear'
  if(maxval(abs(2d0*action_normal_half-action_normal_one))>1d-13) error stop 'normal action is not linear'
  if(maxval(abs(action_value_one(:,:,2,:)))<=maxval(abs(action_value_one(:,:,1,:)))) &
    error stop 'second state was not retained'
  raw_neighbor_value=value
  raw_neighbor_normal=normal
  raw_neighbor_value(:,:,1,:)=-value(:,:,2,:)
  raw_neighbor_value(:,:,2,:)=value(:,:,1,:)
  raw_neighbor_normal(:,:,1,:)=-normal(:,:,2,:)
  raw_neighbor_normal(:,:,2,:)=normal(:,:,1,:)
  projected_neighbor_value=value
  projected_neighbor_normal=normal
  call apply_dg_dc_frozen_projected_face(value,normal,raw_neighbor_value,raw_neighbor_normal,&
    17,17,merge(1,-1,rank==0),0.4d0,0.16d0,9d0,1d0,raw_action_value,raw_action_normal,ok,message)
  if(.not.ok)error stop trim(message)
  call apply_dg_dc_frozen_projected_face(value,normal,projected_neighbor_value,projected_neighbor_normal,&
    17,17,merge(1,-1,rank==0),0.4d0,0.16d0,9d0,1d0,projected_action_value,projected_action_normal,ok,message)
  if(.not.ok)error stop trim(message)
  if(maxval(abs(raw_action_value-projected_action_value))+&
    maxval(abs(raw_action_normal-projected_action_normal))<1d-8)&
    error stop 'same-index transformed action unexpectedly matched projected reference'
  call apply_dg_dc_frozen_projected_face(value,normal,projected_neighbor_value,projected_neighbor_normal,&
    18,17,merge(1,-1,rank==0),0.4d0,0.16d0,9d0,1d0,projected_action_value,projected_action_normal,ok,message)
  if(ok)error stop 'stale frozen projector generation was accepted'
  if(rank==0) print '(a)','PASS direct DC SIPG face MPI fixture'
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif
end program test_dg_dc_direct_cg_mpi
