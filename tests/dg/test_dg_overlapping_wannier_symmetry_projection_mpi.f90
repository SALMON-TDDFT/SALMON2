#include "config.h"
program test_dg_overlapping_wannier_symmetry_projection_mpi
  use mpi
  use dg_overlapping_wannier_symmetry,only:project_dg_fragment_covariant_operators
  implicit none
  integer::ierr,rank,nproc
  complex(8)::representation(2,2,2),scalars(2,2,2),vectors(2,2,3,2)
  real(8)::rotations(3,3,2),pre_defect,post_defect
  complex(8),allocatable::projected_scalars(:,:,:),projected_vectors(:,:,:,:)
  logical::ok
  character(256)::message
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  representation=(0d0,0d0);representation(1,1,1)=1d0;representation(2,2,1)=1d0
  representation(1,1,2)=1d0;representation(2,2,2)=-1d0
  rotations=0d0;rotations(1,1,1)=1d0;rotations(2,2,1)=1d0;rotations(3,3,1)=1d0
  rotations(1,1,2)=-1d0;rotations(2,2,2)=-1d0;rotations(3,3,2)=-1d0
  scalars=(0d0,0d0);scalars(:,:,1)=reshape([cmplx(2d0,0d0,8),cmplx(1d-7,0d0,8), &
    cmplx(1d-7,0d0,8),cmplx(1d0,0d0,8)],[2,2]);scalars(:,:,2)=2d0*scalars(:,:,1)
  vectors=(0d0,0d0)
  vectors(:,:,1,1)=reshape([cmplx(1d-7,0d0,8),cmplx(1d0,0d0,8), &
    cmplx(1d0,0d0,8),cmplx(-1d-7,0d0,8)],[2,2])
  vectors(:,:,2,1)=2d0*vectors(:,:,1,1);vectors(:,:,3,1)=3d0*vectors(:,:,1,1)
  vectors(:,:,:,2)=4d0*vectors(:,:,:,1)
  call project_dg_fragment_covariant_operators(representation,rotations,scalars,vectors,1d-10, &
    projected_scalars,projected_vectors,pre_defect,post_defect,ok,message)
  call require(ok,trim(message));call require(pre_defect>1d-8,'fixture needs nonzero covariance noise')
  call require(post_defect<1d-13,'projected scalar/vector covariance')
  call require(abs(projected_scalars(1,2,1))<1d-14,'scalar inversion projection')
  call require(abs(projected_vectors(1,1,1,1))<1d-14,'polar-vector inversion projection')

  call project_dg_fragment_covariant_operators(representation(:,:,1:1),rotations(:,:,1:1), &
    scalars,vectors,1d-10,projected_scalars,projected_vectors,pre_defect,post_defect,ok,message)
  call require(ok.and.all(projected_scalars==scalars).and.all(projected_vectors==vectors), &
    'C1 projection must be byte preserving')

  scalars(1,2,1)=1d-2;scalars(2,1,1)=1d-2
  call project_dg_fragment_covariant_operators(representation,rotations,scalars,vectors,1d-10, &
    projected_scalars,projected_vectors,pre_defect,post_defect,ok,message)
  call require(.not.ok.and.index(message,'pre-projection')>0,'large covariance defect must fail')
  if(rank==0)write(*,'(a,i0,a)')'PASS overlapping-Wannier symmetry projection on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_failure/=0)error stop label
  end subroutine require
end program test_dg_overlapping_wannier_symmetry_projection_mpi
