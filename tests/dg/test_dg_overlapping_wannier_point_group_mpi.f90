#include "config.h"
program test_dg_overlapping_wannier_point_group_mpi
  use mpi
  use dg_overlapping_wannier_symmetry,only:build_dg_fragment_group_representation
  implicit none
  integer::ierr,rank,nproc
  integer::product_table(2,2)
  complex(8)::metric(2,2),raw(2,2,2)
  complex(8),allocatable::representation(:,:,:)
  real(8)::raw_defect,unitarity_defect,closure_defect
  logical::ok
  character(256)::message
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  product_table=reshape([1,2,2,1],[2,2])
  metric=reshape([cmplx(2d0,0d0,8),cmplx(0.3d0,0d0,8), &
                  cmplx(0.3d0,0d0,8),cmplx(2d0,0d0,8)],[2,2])
  raw=(0d0,0d0);raw(1,1,1)=1.000001d0;raw(2,2,1)=1.000001d0
  raw(1,2,2)=0.999999d0;raw(2,1,2)=0.999999d0
  call build_dg_fragment_group_representation(metric,raw,product_table,1d-10,representation, &
    raw_defect,unitarity_defect,closure_defect,ok,message)
  call require(ok,trim(message))
  call require(raw_defect>1d-7,'raw representation must need polar correction')
  call require(unitarity_defect<1d-12,'metric unitarity after polar correction')
  call require(closure_defect<1d-12,'C2 multiplication closure')
  call require(maxval(abs(representation(:,:,1)-identity2()))<1d-12,'identity representation')

  metric=identity2()
  raw(:,:,1)=identity2()
  raw(:,:,2)=reshape([cmplx(0d0,0d0,8),cmplx(1d0,0d0,8), &
                      cmplx(-1d0,0d0,8),cmplx(0d0,0d0,8)],[2,2])
  call build_dg_fragment_group_representation(metric,raw,product_table,1d-10,representation, &
    raw_defect,unitarity_defect,closure_defect,ok,message)
  call require(.not.ok.and.index(message,'closure')>0,'nonrepresentation must fail closure')
  if(rank==0)write(*,'(a,i0,a)')'PASS overlapping-Wannier point group on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  function identity2()result(identity)
    complex(8)::identity(2,2)
    identity=(0d0,0d0);identity(1,1)=1d0;identity(2,2)=1d0
  end function identity2
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_failure/=0)error stop label
  end subroutine require
end program test_dg_overlapping_wannier_point_group_mpi
