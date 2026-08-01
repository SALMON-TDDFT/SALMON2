#include "config.h"
program test_dg_overlapping_wannier_fragment_symmetry_mpi
  use mpi
  use dg_overlapping_wannier_symmetry, only: select_dg_exact_fragment_subgroup
  implicit none
  integer :: ierr,rank,nproc,i,j
  integer :: product_table(4,4)
  integer :: invalid_product_table(4,4)
  integer,allocatable :: subgroup(:)
  real(8) :: atom(4),boundary(4),grid(4),center(4)
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  do i=1,4
    do j=1,4
      product_table(i,j)=modulo(i+j-2,4)+1
    end do
  end do

  atom=1d-13; boundary=2d-13; grid=3d-13; center=4d-13
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,2,3,4]),'roundoff-scale exact C4')

  atom=[1d-13,2d-4,1d-13,2d-4]
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,3]),'physical displacement reduces C4 to C2')

  atom=[1d-13,2d-4,3d-4,4d-4]
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.size(subgroup)==1.and.subgroup(1)==1,'physical displacement permits C1')

  atom=1d-13; grid(2)=2d-4
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(ok.and.all(subgroup==[1,3]),'grid-incompatible operations form closed subgroup')

  grid=1d-13; boundary(1)=2d-4
  call select_dg_exact_fragment_subgroup(product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(.not.ok.and.index(message,'identity')>0,'invalid identity is fatal')

  invalid_product_table=product_table
  invalid_product_table(2,4)=2
  invalid_product_table(4,2)=2
  boundary=1d-13
  call select_dg_exact_fragment_subgroup(invalid_product_table,atom,boundary,grid,center, &
    1d-10,1d-10,1d-10,1d-10,subgroup,ok,message)
  call require(.not.ok.and.index(message,'group')>0,'non-group product table is fatal')

  if(rank==0)write(*,'(a,i0,a)')'PASS exact buffered-fragment symmetry on ',nproc,' ranks'
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
end program test_dg_overlapping_wannier_fragment_symmetry_mpi
